function HyperLouvain(H::hypergraph,kmax::Int64,Ω,maxits::Int64=100,bigInt::Bool=true;α,verbose=true)
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


    Z = collect(1:n)                    # All nodes start in singleton clusters
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
    V, μ, M = evalSums(Z, H.D, r;constants=false, bigInt=bigInt)
    C = evalConstants(r)

    # Penalties for completely partitioned hypereges
    Pn = zeros(kmax)
    for i = 1:kmax
        p = ones(i)
        om_z = Ω(p; α=α, mode="partition")
        Pn[i] =log(om_z) # later scale this by weight
    end

    edge2penalty = zeros(length(Hyp))    # edge id -> penalty at edge
    for e = 1:length(Hyp)
        edge = Hyp[e]
        weight = w[e]
        k = length(edge)
        edge2penalty[e] = weight*Pn[k]
    end

    # Main Greedy Local Move Loop
    while improving && iter < maxits

        iter += 1
        if mod(iter,1) == 0; 
            if verbose println("Louvain Iteration $iter") end
        end
        improving = false

        for i = 1:n                     # visit each node in turn

            Ci_ind = Z[i]               # Cluster index for node i
            Ci = Clusters[Ci_ind]       # Indices of nodes in i's cluster
            Ni = Neighbs[i]             # Get the indices of i's neighbors
            NC = unique(Z[Ni])          # Clusters adjacent to i (candidate moves)

            BestZ = Z[i]                # The default is do nothing
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
                    ΔV, Δμ, ΔM = increments(V, μ, M, [i], Cj_ind, H.D, Z)
                    voldiff = 0
                    for p in keys(M)
                        voldiff += Ω(p;α=α,mode = "partition")*ΔM[p]*C[p]
                    end
                    cdiff, Δpen = CutDiff_OneNode(Hyp,w,node2edges,Z,i,Cj_ind,edge2penalty,Ω;α=α)
                    change =  cdiff - voldiff
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
                ci_old = Z[i]
                Z[i] = BestZ

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
    Z, Clusters = renumber(Z,Clusters)
    return Z

end

function compute_moddiff(edge2part,Cts,Hyp,w,node2edges,V::Array, μ::Array, M::Dict,I_::T, t::Int64, D::Vector{Int64}, Z::Vector{Int64},C::Dict,Ω;α) where T <: Union{Int64, Vector{Int64}}
    """
    NOTE: I_ can now be a Vector{Int64} of nodes, all of which are assumed to belong in the same cluster.
    """
    # increments due to proposal
    ΔV, Δμ, ΔM = increments(V, μ, M, I_, t, D, Z)

    voldiff = 0
    for p in keys(M)
        voldiff += Ω(p;α=α,mode = "partition")*ΔM[p]*C[p]
    end

    ΔC, Δe2p = CutDiff_Many(Cts, I_, t, Z, Hyp, w, node2edges,edge2part)
    cdiff = sum(ΔC[p]*log(Ω(p; α=α, mode="partition")) for p in keys(ΔC))
    mdiff = cdiff-voldiff

    return mdiff, ΔV, Δμ, ΔM, ΔC, Δe2p
end


function SuperNodeStep(H::hypergraph,Z::Vector{Int64},kmax::Int64,Ω,maxits::Int64=100,bigInt::Bool=true;α,verbose=true)
    """
    A Louvain step, but starting with all nodes in an arbitrary initial cluster
    assignment Z. Louvain only considers moving an entire cluster at once.
    Running this code with Z = collect(1:n) is equivalent to HypergraphLouvain

    H: hypergraph
    Z: initial cluster assingment as a length n vector
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended.
    return: Z::array{Int64, 1}: array of group labels.
    """
    Hyp, w = hyperedge_formatting(H)    # hyperedge to node list
    node2edges = EdgeMap(H)             # node to hyperedge list
    n = length(H.D)

    Z = renumber(Z)     # Ensure clusters go from 1 to number unique clusters
    sn = maximum(Z)     # number of supernodes

    Clusters = Vector{Vector{Int64}}()  # Also store as cluster list
    for j = 1:sn
        push!(Clusters, Vector{Int64}())
    end
    for v = 1:n
        push!(Clusters[Z[v]],v) # put node v in cluster Z[v]
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
    V, μ, M = evalSums(Z, H.D, r;constants=false, bigInt=bigInt);
    C = evalConstants(r)

    # Data structures for fast local changes in cuts
    Cts = evalCuts(Z,H)                    # edge partition -> no. times it appears in Z
    edge2part = Vector{Vector{Int64}}()    # edge id -> partition type
    for e = 1:length(Hyp)
        edge = Hyp[e]
        push!(edge2part,partitionize(Z[edge]))
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
            Ci_ind = Z[rep]         # Cluster associated with this supernode
            Ci = Clusters[Ci_ind]   # Indices of nodes in this cluster

            Ni = Neighbs[i]         # All supernodes that neighbor this supernode

            # Get all the clusters of neighboring supernodes
            NC = Vector{Int64}()
            for ni in Ni
                push!(NC,Z[SuperNodes[ni][1]])
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
                    change, ΔV, Δμ, ΔM, ΔC, Δe2p = compute_moddiff(edge2part,Cts,Hyp, w, node2edges,V, μ, M, S, Cj_ind, H.D, Z, C, Ω; α=α)
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
                ci_old = Z[rep] # index of the old cluster
                Z[S] .= BestZ

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
    Z, Clusters = renumber(Z,Clusters)
    if verbose
        @show t1
    end
    return Z, changemade
end

function SuperNodeLouvain(H::hypergraph,kmax::Int64,Ω,maxits::Int64=100,bigInt::Bool=true;α,verbose=true)
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
    Z = HyperLouvain(H,kmax,Ω;α=α,verbose=verbose)
    n = length(Z)

    changed = false
    if maximum(Z) != n
        changed = true
    end

    while changed
        phase += 1
        if verbose println("SuperNode Louvain: Phase $phase") end
        Z, changed = SuperNodeStep(H,Z,kmax,Ω;α=α,verbose=verbose)
    end

    return Z
end


function HyperLouvain_0(H::hypergraph,kmax::Int64,Ω,maxits::Int64=100,bigInt::Bool=true;α,verbose=verbose)
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


    Z = collect(1:n)                    # All nodes start in singleton clusters
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
    V, μ, M = evalSums(Z, H.D, r;constants=false, bigInt=bigInt)
    C = evalConstants(r)

    # Data structures for fast local changes in cuts
    Cts = evalCuts(Z,H)                    # edge partition -> no. times it appears in Z
    edge2part = Vector{Vector{Int64}}()    # edge id -> partition type
    for e = 1:length(Hyp)
        edge = Hyp[e]
        push!(edge2part,partitionize(Z[edge]))
    end

    # Main Greedy Local Move Loop
    while improving && iter < maxits

        iter += 1
        if mod(iter,1) == 0; 
            if verbose println("Louvain Iteration $iter") end
        end
        improving = false

        for i = 1:n                     # visit each node in turn

            Ci_ind = Z[i]               # Cluster index for node i
            Ci = Clusters[Ci_ind]       # Indices of nodes in i's cluster
            Ni = Neighbs[i]             # Get the indices of i's neighbors
            NC = unique(Z[Ni])          # Clusters adjacent to i (candidate moves)
            # println("$i = i, length(NC) = $(length(NC))")

            BestZ = Z[i]                # The default is do nothing
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
                    change, ΔV, Δμ, ΔM, ΔC, Δe2p = compute_moddiff(edge2part,Cts,Hyp, w, node2edges,V, μ, M, [i], Cj_ind, H.D, Z, C, Ω; α=α)
                end

                # Check if this is currently the best possible greedy move
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

            # Move i to the best new cluster, if it strictly improves modularity
            if BestImprove > 1e-8

                changemade = true
                improving = true    # we have a reason to keep iterating!

                # update increments for volume computation
                V, μ, M = addIncrements(V, μ, M, V_best, μ_best, M_best)

                # update clustering
                ci_old = Z[i]
                Z[i] = BestZ

                # Remove i from its old cluster and add it to its new cluster
                Clusters[ci_old] = setdiff(Clusters[ci_old],i)
                push!(Clusters[BestZ],i)

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

    if ~changemade
        if verbose println("No nodes moved clusters") end
    end
    Z, Clusters = renumber(Z,Clusters)
    return Z

end
