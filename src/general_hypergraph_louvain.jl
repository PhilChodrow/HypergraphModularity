function HyperLouvain(H::hypergraph,kmax::Int64,Ω::IntensityFunction,maxits::Int64=100,bigInt::Bool=true;α,verbose=true,scan_order="lexical", Z0 = collect(1:length(H.D)))
    HyperLouvain(H,Ω;α=α,kmax=kmax,maxits=maxits,bigInt=bigInt,verbose=verbose,scan_order=scan_order, Z0 = Z0)
end

function HyperLouvain(H::hypergraph,Ω::IntensityFunction;α,kmax=maximum(keys(H.E)),maxits::Int64=100,bigInt::Bool=true,verbose=true,scan_order="lexical", Z0 = collect(1:length(H.D)))
    """
    Basic step Louvain algorithm: iterate through nodes and greedily move
    nodes to adjacent clusters. Does not form supernodes and does not recurse.

    H: hypergraph
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended.
    return: Z::array{Int64, 1}: array of group labels.
    """
    Hyp, w = hyperedge_formatting(H)    # hyperedge to node list
    node2edges = EdgeMap(H)             # node to hyperedge list
    n = length(H.D)

    Neighbs = NeighborList(node2edges, Hyp)  # Store neighbors of each node


    Ẑ = Z0                              # All nodes start in singleton clusters
    Clusters = Vector{Vector{Int64}}()  # Also store as cluster list
    for v = 1:n
        push!(Clusters, Vector{Int64}())
        push!(Clusters[v],v)
    end

    iter = 0
    improving = true
    changemade = false

    # Data structures for fast local changes in volumes
    r = kmax
    V, μ, M = evalSums(Ẑ, H.D, r;constants=false, bigInt=bigInt)
    C = evalConstants(r)

    cuts = evalCuts(H, Ẑ, Ω)

    # initialize memoization of edge penalties
    Pn           = [log(Ω.ω(Ω.P(collect(1:i)), α)) for i in 1:kmax]
    edge2penalty = [w[e]*Pn[length(Hyp[e])] for e in 1:length(Hyp)]


    # Main Greedy Local Move Loop
    while improving && iter < maxits

        iter += 1
        if mod(iter,1) == 0;
            if verbose println("Louvain Iteration $iter, Q = $(round(Float64(modularity(H, Ẑ, Ω;α = α)), digits = 2))") end
        end
        improving = false


        scan = scan_order == "lexical" ? (1:n) : Random.shuffle(1:n)

        for i = scan                     # visit each node in turn

            Ci_ind = Ẑ[i]               # Cluster index for node i
            Ci = Clusters[Ci_ind]       # Indices of nodes in i's cluster
            Ni = Neighbs[i]             # Get the indices of i's neighbors
            NC = unique(Ẑ[Ni])          # Clusters adjacent to i (candidate moves)

            BestZ = Ẑ[i]                # The default is do nothing
            BestImprove = 0
            V_best = 0
            μ_best = 0
            M_best = 0
            pen_best = 0

            # Check if it's better to move to a nearby cluster, Cj
            for j = 1:length(NC)
                Cj_ind = NC[j]

                # Check improvement for move from cluster i to to cluster j
                if Cj_ind == Ci_ind
                    change = 0
                else
                    ΔV, Δμ, ΔM = increments(V, μ, M, [i], Cj_ind, H.D, Ẑ)
                    voldiff = sum(Ω.ω(Ω.aggregator(p),α)*ΔM[p]*C[p] for p in keys(ΔM))

                    cdiff, Δpen = CutDiff_OneNode(Hyp,w,node2edges,Ẑ,i,Cj_ind,edge2penalty,Ω;α=α) # NOTE: need to check on this
                    change =  cdiff - voldiff
#                     println(cdiff, " ", voldiff)
                end

                # Check if this is currently the best possible greedy move
                if change - BestImprove > 1e-8
                    V_best = ΔV
                    μ_best = Δμ
                    M_best = ΔM
                    pen_best = Δpen
                    BestImprove = change
                    BestZ = Cj_ind
                    improving = true
                end
            end

            # Move i to the best new cluster, if it strictly improves modularity
            if BestImprove > 1e-8

                changemade = true
                improving = true    # we have a reason to keep iterating!

                # update increments for volume computation
                V, μ, M = addIncrements(V, μ, M, V_best, μ_best, M_best)

                # update clustering
                ci_old = Ẑ[i]
                Ẑ[i] = BestZ

                # Remove i from its old cluster and add it to its new cluster
                Clusters[ci_old] = setdiff(Clusters[ci_old],i)
                push!(Clusters[BestZ],i)

                # update map from edge ID to penalty at edge
                for eid in keys(pen_best)
                    newpenalty = pen_best[eid]
                    edge2penalty[eid] = newpenalty
                end
            end
        end
    end

    if ~changemade
        if verbose println("No nodes moved clusters") end
    end
    Ẑ, Clusters = renumber(Ẑ,Clusters)
    return Ẑ

end

function compute_moddiff(edge2part,Cts,Hyp,w,node2edges,V::Array, μ::Array, M::Dict,I_::T, t::Int64, D::Vector{Int64}, Z::Vector{Int64},C::Dict,Ω::IntensityFunction;α) where T <: Union{Int64, Vector{Int64}}
    """
    NOTE: I_ can now be a Vector{Int64} of nodes, all of which are assumed to belong in the same cluster.
    """
    # increments due to proposal
    ΔV, Δμ, ΔM = increments(V, μ, M, I_, t, D, Z)

    voldiff = sum(Ω.ω(Ω.aggregator(p),α)*ΔM[p]*C[p] for p in keys(M))

    ΔC, Δe2p = CutDiff_Many(Cts, I_, t, Z, Hyp, w, node2edges,edge2part)
    cdiff = sum(ΔC[p]*log(Ω.ω(Ω.aggregator(p), α)) for p in keys(ΔC))
    mdiff = cdiff-voldiff

    return mdiff, ΔV, Δμ, ΔM, ΔC, Δe2p
end


function SuperNodeStep(H::hypergraph,Z::Vector{Int64},kmax::Int64,Ω,maxits::Int64=100,bigInt::Bool=true;α,verbose=true)
    """
    A Louvain step, but starting with all nodes in an arbitrary initial cluster
    assignment Z. Louvain only considers moving an entire cluster at once.
    Running this code with Z = collect(1:n) is equivalent to HypergraphLouvain

    H: hypergraph
    Z: initial cluster assignment as a length n vector
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended.
    return: Z::array{Int64, 1}: array of group labels.
    """
    Hyp, w = hyperedge_formatting(H)    # hyperedge to node list
    node2edges = EdgeMap(H)             # node to hyperedge list
    n = length(H.D)

    Ẑ = renumber(Z)     # Ensure clusters go from 1 to number unique clusters
    sn = maximum(Ẑ)     # number of supernodes

    Clusters = Vector{Vector{Int64}}()  # Also store as cluster list
    for j = 1:sn
        push!(Clusters, Vector{Int64}())
    end
    for v = 1:n
        push!(Clusters[Ẑ[v]],v) # put node v in cluster Z[v]
    end
    SuperNodes = deepcopy(Clusters)

    # Store indices of supernodes that neighbor each supernode
    tic = time()
    Neighbs = SuperNeighborList(Hyp, SuperNodes,n)
    toc = time()-tic
    if verbose
        @show toc
    end

    # Note: There's a very subtle difference between the initial clustering,
    # which defines a set of supernodes, and the current clustering that
    # is being built by Louvain.
    # Supernodes don't change--they move as a unit. Clusters change.

    iter = 0
    improving = true
    changemade = false

    # Data structures for fast local changes in volumes
    r = kmax
    V, μ, M = evalSums(Ẑ, H.D, r;constants=false, bigInt=bigInt);
    C = evalConstants(r)

    # Data structures for fast local changes in cuts
    Cts = evalCuts(Ẑ,H)                    # edge partition -> no. times it appears in Z
    edge2part = Vector{Vector{Int64}}()    # edge id -> partition type
    for e = 1:length(Hyp)
        edge = Hyp[e]
        push!(edge2part,partitionize(Ẑ[edge]))
    end
    t1 = 0
    # Main Greedy Local Move Loop
    while improving && iter < maxits

        iter += 1
        if mod(iter,1) == 0;
            if verbose println("Louvain Iteration $iter") end
        end
        improving = false


        for i = 1:sn          # visit each supernode in turn

            S = SuperNodes[i]       # Indices associated with this supernode
            rep = S[1]              # A "representative" node from this supernode
            Ci_ind = Ẑ[rep]         # Cluster associated with this supernode
            Ci = Clusters[Ci_ind]   # Indices of nodes in this cluster

            Ni = Neighbs[i]         # All supernodes that neighbor this supernode

            # Get all the clusters of neighboring supernodes
            NC = Vector{Int64}()
            for ni in Ni
                push!(NC,Ẑ[SuperNodes[ni][1]])
            end
            NC = unique(NC)
            # println("$i = i, length(NC) = $(length(NC))")

            BestZ = Ci_ind          # The default is to not move the cluster
            BestImprove = 0
            V_best = 0
            μ_best = 0
            M_best = 0
            C_best = 0
            e2p_best = 0

            # Check if it's better to move to a nearby cluster, Cj
            for j = 1:length(NC)
                Cj_ind = NC[j]

                # Check improvement for move from cluster i to to cluster j
                if Cj_ind == Ci_ind
                    change = 0
                else
                    tic = time()
                    change, ΔV, Δμ, ΔM, ΔC, Δe2p = compute_moddiff(edge2part,Cts,Hyp, w, node2edges,V, μ, M, S, Cj_ind, H.D, Ẑ, C, Ω; α=α)
                    t1 += time()-tic
                end

                # Check if this is currently the best possible greedy move to make
                if change - BestImprove > 1e-8
                    V_best = ΔV
                    μ_best = Δμ
                    M_best = ΔM
                    C_best = ΔC
                    e2p_best = Δe2p
                    BestImprove = change
                    BestZ = Cj_ind
                    improving = true
                end
            end

            # Move supernode i to the best new cluster, only if it strictly improves modularity
            if BestImprove > 1e-8

                changemade = true
                improving = true # we have a reason to keep iterating!

                # update increments for volume computation
                V, μ, M = addIncrements(V, μ, M, V_best, μ_best, M_best)

                # update clustering
                ci_old = Ẑ[rep] # index of the old cluster
                Ẑ[S] .= BestZ

                # Remove all nodes from supernode i from their old cluster.
                Clusters[ci_old] = setdiff(Clusters[ci_old],S)

                # Add them to their new cluster
                append!(Clusters[BestZ],S)

                # Update map from partition type to number of those partitions
                for part in keys(C_best)
                    Δcount = C_best[part]
                    Cts[part] = get(Cts,part,0) + Δcount
                end

                # update map from edge ID to partition type
                for eid in keys(e2p_best)
                    newpart = e2p_best[eid]
                    edge2part[eid] = newpart
                end
            end
        end
    end
    Ẑ, Clusters = renumber(Ẑ,Clusters)
    if verbose
        @show t1
    end
    return Ẑ, changemade
end

function SuperNodeLouvain(H::hypergraph,kmax::Int64,Ω,maxits::Int64=100,bigInt::Bool=true;α,verbose=true,scan_order="lexical", Z0 = nothing)
    SuperNodeLouvain(H,Ω;α=α,kmax=kmax,maxits=maxits,bigInt=bigInt,verbose=verbose,scan_order=scan_order, Z0 = copy(Z0))
end

function SuperNodeLouvain(H::hypergraph,Ω::IntensityFunction;α,kmax=maximum(keys(H.E)),maxits::Int64=100,bigInt::Bool=true,verbose=true,scan_order="lexical", Z0 = collect(1:length(H.D)))
    """
    Running Louvain and then the super-node louvain steps until no more
    progress is possible

    H: hypergraph
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended.
    return: Z::array{Int64, 1}: array of group labels.
    """

    phase = 1
    
    if verbose println("Faster SuperNode Louvain: Phase $phase") end
    Ẑ = HyperLouvain(H,kmax,Ω;α=α,verbose=verbose,scan_order=scan_order, Z0 = copy(Z0))
    
    
    n = length(Ẑ)

    changed = true

    while changed
        phase += 1
        if verbose println("SuperNode Louvain: Phase $phase") end
        Ẑ, changed = SuperNodeStep(H,Ẑ,kmax,Ω,maxits;α=α,verbose=verbose)
    end

    return Ẑ
end


function evalCuts(Z::Array{Int64,1}, H::hypergraph)
    """
    Gets map from partition type to the number of times that partition type
    shows up in edges for the clustering Z
    """
    C = Dict{Vector{Int64}, Int64}()
    for k in keys(H.E)
        Ek = H.E[k]
        for e in keys(Ek)
            p = partitionize(Z[e])
            C[p] = get(C, p, 0) + Ek[e]
        end
    end
    return C
end

function evalCuts(H::hypergraph, Z::Array{Int64,1}, Ω::IntensityFunction)

    C = Dict(p => 0 for p in Ω.range)
    for k in keys(H.E)
        Ek = H.E[k]
        for e in keys(Ek)
            p = Ω.P(Z[e])
            C[p] = get(C, p, 0) + Ek[e]
        end
    end
    return C
end

function CutDiff_Many(C, I, t, Z, Hyp, w, node2edges,edge2part)
    """
    C: Dict{Vector{Int64}, Int64}, the dict of current cut values as would be produced by evalCuts().
        C is keyed by partition vectors p.
        C[p] is the # of edges with label partition p.
        Note: Ω is NOT used in the calculation of C
    I: Vector{Int64}, the list of nodes to move, assumed to be in the same group
    t: Int64, the new group to which to move the nodes I
    Z: Vector{Int64}, the current group assignments
    Hyp: list of nodes in each hyperedge
    w:  weight of the hyperedge (number of times this set of nodes is assigned a hyperedge)
    node2edges: node ID to list of edges the node is in.
    edge2part: for hyperedge index e, returns the current partitionize vector for e for the clustering Z
    return: ΔC, a Dict{Vector{Int64}, Int64}, where ΔC[p] is the change in C[p] caused by moving the nodes in I to group t
    """

    E_id = Vector{Int64}()  # edge IDs that invole some node from I
    for i in I
        append!(E_id,node2edges[i])
    end
    E = unique(E_id)

    # This stores proposed changes to edge partitions
    Δe2p = Dict{Int64,Vector{Int64}}()

    # Store changes in the number of partition types in the clustering
    ΔC = Dict(p => 0 for p in keys(C))

    Z_prop = copy(Z)    # Store the proposed new clustering
    Z_prop[I] .= t

    for eid in E
        e = Hyp[eid]    # get nodes in edge
        we = w[eid]     # get the edge weight

        p_old = edge2part[eid]  # instead of: pold = partitionize(Z[e])
        p_new = partitionize(Z_prop[e])
        Δe2p[eid] = p_new

        # Update number of partition types for p_old and p_new
        ΔC[p_old] -= we
        ΔC[p_new] = get(ΔC,p_new,0) + we

    end
    return ΔC, Δe2p
end


function CutDiff_OneNode(H::Vector{Vector{Int64}},w::Array{Float64,1},node2edges::Vector{Vector{Int64}},Z::Array{Int64,1}, I::Int64,J::Int64,edge2penalty::Union{Vector{Float64}, Vector{BigFloat}},Ω::IntensityFunction; α)
    """
    CutDiff: Compute change in the first term of the modularity function
    resulting from moving a node I to cluster J.
    """
    orig = Z[I]
    He   = node2edges[I]
    Δpen = Dict{Int64,Float64}()    # Change in penalty for moving a node

    obj  = 0
    Z[I] = J
    for i in He
        e        = H[i]
        z        = Z[e]                      # set of clusters
        new_pen = w[i]*log(Ω.ω(Ω.P(z),α))   # new penalty if this move is made
        Δpen[i]  = new_pen
        obj     += Δpen[i]-edge2penalty[i]
    end
    Z[I] = orig

    return obj, Δpen
end

function SuperNeighborList(Hyp, SuperNodes,n)
    """
    For an initial clustering into supernodes, get a list of other supernodes
    adjacenct to each supernode
    """
    H = elist2incidence(Hyp,n)
    m = size(H,1)
    sn = length(SuperNodes)
    I = Vector{Int64}()
    J = Vector{Int64}()

    for s = 1:sn
        # For supernode j, get indices of all other supernodes i such that
        # some node in j shares a node with some node in i
        S = SuperNodes[s]

        # this is silly and slow. But also not a bottleneck, so I'm
        # not worrying about making it faster right now
        Edges = findall(x->x>0,vec(sum(H[:,S],dims = 2)))  # all hyperedges nodes from S touch

        for e in Edges
            push!(I,s)
            push!(J,e)
        end
    end

    # binary (unweighted) supernodes-by-edges incidence matrix
    Hnew = SparseArrays.sparse(J,I,ones(length(I)),m,sn)
    Neighbs = NeighborList(Hnew, SparseArrays.sparse(Hnew'))
    return Neighbs
end
