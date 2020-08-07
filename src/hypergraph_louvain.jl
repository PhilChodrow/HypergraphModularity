using Random
using SparseArrays

include("objectives.jl")
include("hyper_format.jl")
include("HSBM.jl")
include("diffs.jl")

function cutdiff(He2n::SparseMatrixCSC{Float64,Int64},Hn2e::SparseMatrixCSC{Float64,Int64},w::Array{Float64,1},c::Array{Int64,1}, I::Int64,J::Int64,Ω;α)
    """
    CutDiff: Compute change in the first term of the modularity function
    resulting from moving a node I to cluster J.

    This code is sloppy and I still think there's a faster way, but this is a step.
    """
    orig = c[I]
    He = getedges(He2n,I)
    obj1 = 0
    fprintf("using this code but shouldn't be")
    for i = He
        edge = getnodes(Hn2e,i)
        clus_e = c[edge]    # set of clusters
        p = partitionize(clus_e)
        om_z = Ω(p; α=α, mode="partition")
        obj1 += w[i]*log(om_z)
    end

    c[I] = J
    obj2 = 0
    for i = He
        edge = getnodes(Hn2e,i)
        clus_e = c[edge]    # set of clusters
        p = partitionize(clus_e)
        om_z = Ω(p; α=α,mode="partition")
        obj2 += w[i]*log(om_z)
    end
    c[I] = orig
    return obj2-obj1
end

function renumber(Z::Vector{Int64},Clusters::Vector{Vector{Int}})
    """
    Renumber a cluster vector so that it goes from 1 to # of clusters

    e.g. Z = [2 3 5 3 9] ---> Z = [1 2 3 2 4]

    Z[i] = (integer) cluster ID that node i belongs to
    Clusters[j] = (integer array) node IDs for nodes in cluster j
    """
    n = length(Z)
    map = unique(Z)
    cnew = zeros(Int64,n)

    Clusters = Clusters[map]

    # Rename the clusters
    for i = 1:n
        newClus = findfirst(x->x == Z[i],map)
        cnew[i] = newClus
        push!(Clusters[newClus],i)
    end

    return cnew, Clusters

end


function renumber(c::Vector{Int64})
    """
    See above function, this function doesn't care about the Clusters array
    """
    n = length(c)
    map = sort(unique(c))
    cnew = zeros(Int64,n)

    # Rename the clusters
    for i = 1:n
        newClus = findfirst(x->x == c[i],map)
        cnew[i] = newClus
    end

    return cnew

end

function Naive_HyperLouvain(H::hypergraph,Ω,maxits::Int64=100,bigInt::Bool=true;α)
    """
    Basic step Louvain algorithm: iterate through nodes and greedily move
    nodes to adjacent clusters. Does not form supernodes and does not recurse.

    H: hypergraph
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended.
    return: Z::array{Int64, 1}: array of group labels.
    """
    println("")
    # Begin by converting to alternative hypergraph storage
    He2n, w = hypergraph2incidence(H)
    Hn2e = sparse(He2n')
    m, n = size(He2n)

    # Store node neighbors of each node
    Neighbs = NeighborList(He2n,Hn2e)

    # All nodes start in singleton clusters
    Z = collect(1:n)

    # Also store as cluster list
    Clusters = Vector{Vector{Int64}}()
    for v = 1:n
        push!(Clusters, Vector{Int64}())
        push!(Clusters[v],v)
    end

    improving = true
    iter = 0
    changemade = false

    while improving && iter < maxits

        iter += 1
        if mod(iter,1) == 0
            println("Louvain Iteration $iter")
        end
        improving = false

        # visit each node in turn
        for i = 1:n

            # Cluster index for node i
            Ci_ind = Z[i]

            # Get the indices of nodes in i's cluster
            Ci = Clusters[Ci_ind]

            # Get the indices of i's neighbors--these define clusters we might move to
            Ni = Neighbs[i]
            # Get the neighboring clusters of i:
            # there are the only places we would even consider moving node i to
            NC = unique(Z[Ni])

            # The default is to not move the cluster
            BestC = Z[i]
            BestImprove = 0

            # Now let's see if it's better to move to a nearby cluster, Cj
            for j = 1:length(NC)
                Cj_ind = NC[j]

                # Check how much it would improve to move i to to cluster j
                if Cj_ind == Ci_ind
                    change = 0
                else

                    Znew = copy(Z)
                    Znew[i] = Cj_ind
                    change = modularity(H,Znew,Ω;α=α) - modularity(H,Z,Ω;α=α)

                end

                # Check if this is currently the best possible greedy move to make
                if change - BestImprove > 0
                    BestImprove = change
                    BestC = Cj_ind
                    improving = true
                end
            end

            # Move i to the best new cluster, only if it strictly improves modularity
            if BestImprove > 0

                ci_old = Z[i]
                Z[i] = BestC

                # Remove i from its old cluster...
                Clusters[ci_old] = setdiff(Clusters[ci_old],i)

                # ...and add it to its new cluster
                push!(Clusters[BestC],i)
                changemade = true
                improving = true # we have a reason to keep iterating!
            end
        end
    end

    if ~changemade
        println("No nodes moved clusters")
    end
    Z, Clusters = renumber(Z,Clusters)
    return Z

end


function HyperLouvain(H::hypergraph,kmax::Int64,Ω,maxits::Int64=100,bigInt::Bool=true;α, verbose=true)
    """
    Basic step Louvain algorithm: iterate through nodes and greedily move
    nodes to adjacent clusters. Does not form supernodes and does not recurse.

    H: hypergraph
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended.
    return: Z::array{Int64, 1}: array of group labels.
    """
    # Begin by converting to alternative hypergraph storage
    start = time()
    Hyp, w = hyperedge_formatting(H)    # hyperedge to node list
    node2edges = EdgeMap(H)             # node to hyperedge list
    convert = time()-start
    n = length(H.D)
    if verbose println("") end
    # println("Took $convert seconds to convert to different hypergraph formats")

    # Store node neighbors of each node
    Neighbs = NeighborList(node2edges, Hyp)

    # All nodes start in singleton clusters
    Z = collect(1:n)

    # Also store as cluster list
    Clusters = Vector{Vector{Int64}}()
    for v = 1:n
        push!(Clusters, Vector{Int64}())
        push!(Clusters[v],v)
    end

    improving = true
    iter = 0
    changemade = false

    # Code for fast local changes in volumes
    r = kmax
    V, μ, M = evalSums(Z, H.D, r;constants=false, bigInt=bigInt);
    C = evalConstants(r)
    while improving && iter < maxits

        iter += 1
        if mod(iter,1) == 0
            if verbose println("Louvain Iteration $iter") end
        end
        improving = false

        # visit each node in turn
        for i = 1:n

            # println("NODE $i")

            # Cluster index for node i
            Ci_ind = Z[i]

            # Get the indices of nodes in i's cluster
            Ci = Clusters[Ci_ind]

            # Get the indices of i's neighbors--these define clusters we might move to
            Ni = Neighbs[i]
            # Get the neighboring clusters of i:
            # there are the only places we would even consider moving node i to
            NC = unique(Z[Ni])

            # The default is to not move the cluster
            BestC = Z[i]
            BestImprove = 0
            V_best = 0
            μ_best = 0
            M_best = 0

            # Now let's see if it's better to move to a nearby cluster, Cj
            for j = 1:length(NC)
                Cj_ind = NC[j]
                # Check how much it would improve to move i to to cluster j
                if Cj_ind == Ci_ind
                    change = 0
                else
                    # Naive evaluation for debugging...
                    # Znew = copy(Z)
                    # Znew[i] = Cj_ind
                    # voldiff1 = second_term_eval(H, Znew, Ω; bigInt = bigInt,α=α )-second_term_eval(H, Z, Ω; bigInt = bigInt, α=α)
                    # cdiff1 = NaiveCutDiff(H,Z,i,Cj_ind, Ω; α=α)

                    voldiff, ΔV, Δμ, ΔM = compute_voldiff2(V, μ, M, i, Cj_ind, H.D, Z, C, Ω; α=α)

                    cdiff = CutDiff(Hyp,w,node2edges,Z,i,Cj_ind,Ω;α=α)

                    change =  cdiff - voldiff
                end

                # Check if this is currently the best possible greedy move to make
                # The tolerance helps with making some behavior more stable:
                #   you only move if there's a numerically nonzero reason to move.
                if change - BestImprove > 1e-8
                    V_best = ΔV
                    μ_best = Δμ
                    M_best = ΔM
                    BestImprove = change
                    BestC = Cj_ind
                    improving = true
                end
            end

            # Move i to the best new cluster, only if it strictly improves modularity
            if BestImprove > 1e-8
                # increments
                V, μ, M = addIncrements(V, μ, M, V_best, μ_best, M_best)

                # update clustering
                ci_old = Z[i]
                Z[i] = BestC

                # Remove i from its old cluster...
                Clusters[ci_old] = setdiff(Clusters[ci_old],i)

                # ...and add it to its new cluster
                push!(Clusters[BestC],i)
                changemade = true
                improving = true # we have a reason to keep iterating!
            end
        end
    end

    if ~changemade
        println("No nodes moved clusters")
    end
    Z, Clusters = renumber(Z,Clusters)
    return Z

end


function compute_voldiff(V::Array, μ::Array, M::Dict,I_::T, t::Int64, D::Vector{Int64}, Z::Vector{Int64},C::Dict,Ω;α) where T <: Union{Int64, Vector{Int64}}
    """
    NOTE: I_ can now be a Vector{Int64} of nodes, all of which are assumed to belong in the same cluster.
    I have not tested compute_voldiff, but I have tested increments() with this usage.
    """
    # increments due to proposal
    ΔV, Δμ, ΔM = increments(V, μ, M, I_, t, D, Z)

    # new proposed quantities
    V_prop, μ_prop, M_prop = addIncrements(V, μ, M, ΔV, Δμ, ΔM)

    vol = 0
    for p in keys(M)
        vol += Ω(p;α=α,mode = "partition")*M[p]*C[p]
    end
    vol_prop = 0
    for p in keys(M_prop)
        vol_prop += Ω(p;α=α,mode = "partition")*M_prop[p]*C[p]
    end
    voldiff = vol-vol_prop

    return voldiff, ΔV, Δμ, ΔM
end

function compute_voldiff2(V::Array, μ::Array, M::Dict,I_::T, t::Int64, D::Vector{Int64}, Z::Vector{Int64},C::Dict,Ω;α) where T <: Union{Int64, Vector{Int64}}
    """
    NOTE: I_ can now be a Vector{Int64} of nodes, all of which are assumed to belong in the same cluster.
    I have not tested compute_voldiff, but I have tested increments() with this usage.
    should give identical behavior as compute_voldiff, has not yet been tested.
    """
    # increments due to proposal
    ΔV, Δμ, ΔM = increments(V, μ, M, I_, t, D, Z)

    # new proposed quantities
    # V_prop, μ_prop, M_prop = addIncrements(V, μ, M, ΔV, Δμ, ΔM)

    voldiff = 0
    for p in keys(M)
        voldiff += Ω(p;α=α,mode = "partition")*ΔM[p]*C[p]
    end

    return voldiff, ΔV, Δμ, ΔM
end



function compute_moddiff(Cts,Hyp,w,node2edges,V::Array, μ::Array, M::Dict,I_::T, t::Int64, D::Vector{Int64}, Z::Vector{Int64},C::Dict,Ω;α) where T <: Union{Int64, Vector{Int64}}
    """
    NOTE: I_ can now be a Vector{Int64} of nodes, all of which are assumed to belong in the same cluster.
    """
    # increments due to proposal
    ΔV, Δμ, ΔM = increments(V, μ, M, I_, t, D, Z)

    voldiff = 0
    for p in keys(M)
        voldiff += Ω(p;α=α,mode = "partition")*ΔM[p]*C[p]
    end

    # are the keys of M and the keys of ΔC always the same? If so, we could
    # maybe combine these pieces to make it faster

    ΔC = cutDiff(Cts, I_, t, Z, Hyp, w, node2edges)
    cdiff = sum(ΔC[p]*log(Ω(p; α=α, mode="partition")) for p in keys(ΔC))
    mdiff = cdiff-voldiff

    return mdiff, ΔV, Δμ, ΔM, ΔC
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
    Hnew = sparse(J,I,ones(length(I)),m,sn)
    Neighbs = NeighborList(Hnew, sparse(Hnew'))
    return Neighbs
end

function Naive_SuperNodeStep(H::hypergraph,Z::Vector{Int64},kmax::Int64,Ω,maxits::Int64=100,bigInt::Bool=true;α)
    """
    A Louvain step, but starting with all nodes in an arbitrary initial cluster
    assignment Z. Louvaint only considers moving an entire cluster at once.
    Running this code with Z = collect(1:n) is equivalent to HypergraphLouvain

    This is a naive first version of the code, which doesn't do anything
    smart about being careful which clusters to even consider moving to.

    H: hypergraph
    Z: initial cluster assingment as a length n vector
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended.
    return: Z::array{Int64, 1}: array of group labels.
    """
    # Begin by converting to alternative hypergraph storage
    Hyp, w = hyperedge_formatting(H)    # hyperedge to node list
    node2edges = EdgeMap(H)             # node to hyperedge list
    n = length(H.D)
    println("")

    # All nodes start in singleton clusters
    Z = renumber(Z)

    sn = maximum(Z) # number of supernodes

    # Also store as cluster list
    Clusters = Vector{Vector{Int64}}()
    for j = 1:sn
        push!(Clusters, Vector{Int64}())
    end
    for v = 1:n
        push!(Clusters[Z[v]],v) # put node v in cluster Z[v]
    end
    SuperNodes = deepcopy(Clusters)

    # Store indices of supernodes that neighbor each supernode
    Neighbs = SuperNeighborList(Hyp, SuperNodes,n)

    # There's a very subtle difference between the initial clustering,
    # which defines a set of supernodes, and the current clustering that
    # is being built by Louvain. Supernodes don't change--they move as a unit.
    # Clusters change.

    improving = true
    iter = 0
    changemade = false

    while improving && iter < maxits

        iter += 1
        if mod(iter,1) == 0
            println("Louvain Iteration $iter")
        end
        improving = false

        # visit each supernode in turn
        for i = 1:sn

            # Get indices associated with this supernode
            S = SuperNodes[i]

            # Grab a "representative" node from this supernode
            rep = S[1]

            # Get the cluster associated with this supernode
            Ci_ind = Z[rep]

            # Get the indices of nodes in this cluster
            Ci = Clusters[Ci_ind]

            # Get all supernodes that neighbor this supernode
            Ni = Neighbs[i]

            # Then get all the clusters of neighboring supernodes
            NC = Vector{Int64}()
            for ni in Ni
                push!(NC,Z[SuperNodes[ni][1]])
            end

            # do something naive first: just check move to all other clusters
            # NC = collect(1:length(Clusters))

            # The default is to not move the cluster
            BestC = Ci_ind
            BestImprove = 0

            # Now let's see if it's better to move to a nearby cluster, Cj
            for j = 1:length(NC)
                Cj_ind = NC[j]

                # Check how much it would improve to move i to to cluster j
                if Cj_ind == Ci_ind
                    change = 0
                else


                    Znew = copy(Z)
                    # Move the entire supernode to the new cluster
                    Znew[S] .= Cj_ind

                    # PC -- look for performance improvements here.
                    change = modularity(H,Znew,Ω;α=α) - modularity(H,Z,Ω;α=α)

                end

                # Check if this is currently the best possible greedy move to make
                if change > BestImprove
                    BestImprove = change
                    BestC = Cj_ind
                    improving = true
                end
            end

            # Move supernode i to the best new cluster, only if it strictly improves modularity
            if BestImprove > 0

                # update clustering
                ci_old = Z[rep] # index of the old cluster
                Z[S] .= BestC

                # Remove all nodes from supernode i from their old cluster.

                Clusters[ci_old] = setdiff(Clusters[ci_old],S)

                # Add them to their new cluster
                append!(Clusters[BestC],S)
                changemade = true
                improving = true # we have a reason to keep iterating!
            end
        end
    end
    Z, Clusters = renumber(Z,Clusters)
    return Z, changemade
end

function SuperNodeStep(H::hypergraph,Z::Vector{Int64},kmax::Int64,Ω,maxits::Int64=100,bigInt::Bool=true;α,verbose=true)
    """
    A Louvain step, but starting with all nodes in an arbitrary initial cluster
    assignment Z. Louvain only considers moving an entire cluster at once.
    Running this code with Z = collect(1:n) is equivalent to HypergraphLouvain

    This is a naive first version of the code, which doesn't do anything
    smart about being careful which clusters to even consider moving to.

    H: hypergraph
    Z: initial cluster assingment as a length n vector
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended.
    return: Z::array{Int64, 1}: array of group labels.
    """
    # Begin by converting to alternative hypergraph storage
    Hyp, w = hyperedge_formatting(H)    # hyperedge to node list
    node2edges = EdgeMap(H)             # node to hyperedge list
    n = length(H.D)
    if verbose println("") end

    # All nodes start in singleton clusters
    Z = renumber(Z)

    sn = maximum(Z) # number of supernodes

    # Also store as cluster list
    Clusters = Vector{Vector{Int64}}()
    for j = 1:sn
        push!(Clusters, Vector{Int64}())
    end
    for v = 1:n
        push!(Clusters[Z[v]],v) # put node v in cluster Z[v]
    end
    SuperNodes = deepcopy(Clusters)

    # Store indices of supernodes that neighbor each supernode
    Neighbs = SuperNeighborList(Hyp, SuperNodes,n)

    # There's a very subtle difference between the initial clustering,
    # which defines a set of supernodes, and the current clustering that
    # is being built by Louvain. Supernodes don't change--they move as a unit.
    # Clusters change.

    improving = true
    iter = 0
    changemade = false

    # Code for fast local changes in volumes
    r = kmax
    V, μ, M = evalSums(Z, H.D, r;constants=false, bigInt=bigInt);
    C = evalConstants(r)
    Cts = evalCuts(Z,H)

    while improving && iter < maxits

        iter += 1
        if mod(iter,1) == 0
            if verbose println("Louvain Iteration $iter") end
        end
        improving = false

        # visit each supernode in turn
        for i = 1:sn

            # Get indices associated with this supernode
            S = SuperNodes[i]

            # Grab a "representative" node from this supernode
            rep = S[1]

            # Get the cluster associated with this supernode
            Ci_ind = Z[rep]

            # Get the indices of nodes in this cluster
            Ci = Clusters[Ci_ind]

            # Get all supernodes that neighbor this supernode
            Ni = Neighbs[i]

            # Then get all the clusters of neighboring supernodes
            NC = Vector{Int64}()
            for ni in Ni
                push!(NC,Z[SuperNodes[ni][1]])
            end

            # The default is to not move the cluster
            BestC = Ci_ind
            BestImprove = 0
            V_best = 0
            μ_best = 0
            M_best = 0
            C_best = 0

            # Now let's see if it's better to move to a nearby cluster, Cj
            for j = 1:length(NC)
                Cj_ind = NC[j]

                # Check how much it would improve to move i to to cluster j
                if Cj_ind == Ci_ind
                    change = 0
                else
                    change, ΔV, Δμ, ΔM, ΔC = compute_moddiff(Cts,Hyp, w, node2edges,V, μ, M, S, Cj_ind, H.D, Z, C, Ω; α=α)
                end

                # Check if this is currently the best possible greedy move to make
                if change > BestImprove
                    V_best = ΔV
                    μ_best = Δμ
                    M_best = ΔM
                    C_best = ΔC
                    BestImprove = change
                    BestC = Cj_ind
                    improving = true
                end
            end

            # Move supernode i to the best new cluster, only if it strictly improves modularity
            if BestImprove > 0

                # increments
                V, μ, M = addIncrements(V, μ, M, V_best, μ_best, M_best)

                # update clustering
                ci_old = Z[rep] # index of the old cluster
                Z[S] .= BestC

                # Remove all nodes from supernode i from their old cluster.
                Clusters[ci_old] = setdiff(Clusters[ci_old],S)

                # Add them to their new cluster
                append!(Clusters[BestC],S)
                changemade = true
                improving = true # we have a reason to keep iterating!

                # Naive update of cuts
                # Phil-- is there an analog of "addIncrements" but for cuts?
                Cts = evalCuts(Z,H)
            end
        end
    end
    Z, Clusters = renumber(Z,Clusters)
    return Z, changemade
end

function SuperNodeLouvain(H::hypergraph,kmax::Int64,Ω,maxits::Int64=100,bigInt::Bool=true;α,verbose=True)
    """
    Running Louvain and then the super-node louvain steps until no more
    progress is possible

    H: hypergraph
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended.
    return: Z::array{Int64, 1}: array of group labels.
    """

    phase = 1
    if verbose println("SuperNode Louvain: Phase $phase") end
    Z = HyperLouvain(H,kmax,Ω;α=α, verbose=verbose)
    n = length(Z)

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

function naiveModDiff(S, j, H, Z, Ω; α)
    """
    a naive computation of the change in modularity associated with changing all the nodes in S from their current cluster (assumed to be the same) to cluster j.
    Simply computes the modularities of the old and new partitions.
    Should be identical to what's currently used at line 487.
    """
    Znew = copy(Z)
    Znew[S] .= Cj_ind
    change = modularity(H,Znew,Ω;α=α) - modularity(H,Z,Ω;α=α)
    return change
end
