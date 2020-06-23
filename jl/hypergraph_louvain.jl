using Random
using SparseArrays

include("objectives.jl")
include("hyper_format.jl")
include("HSBM.jl")

function cutdiff(He2n::SparseMatrixCSC{Float64,Int64},Hn2e::SparseMatrixCSC{Float64,Int64},w::Array{Float64,1},c::Array{Int64,1}, I::Int64,J::Int64,Ω)
    """
    CutDiff: Compute change in the first term of the modularity function
    resulting from moving a node I to cluster J.

    This code is sloppy and I still think there's a faster way, but this is a step.
    """
    orig = c[I]
    He = getedges(He2n,I)
    obj1 = 0
    for i = He
        edge = getnodes(Hn2e,i)
        l = length(edge)
        clus_e = c[edge]    # set of clusters
        p = partitionize(clus_e)
        om_z = Ω(p; mode="partition", k = l)
        obj1 += w[i]*log(om_z)
    end

    c[I] = J
    obj2 = 0
    for i = He
        edge = getnodes(Hn2e,i)
        l = length(edge)
        clus_e = c[edge]    # set of clusters
        p = partitionize(clus_e)
        om_z = Ω(p; mode="partition", k = l)
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

function Naive_HyperLouvain(H::hypergraph,Ω,maxits::Int64=100,bigInt::Bool=true)
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

            # Now let's see if it's better to move to a nearby cluster, Cj
            for j = 1:length(NC)
                Cj_ind = NC[j]

                # Check how much it would improve to move i to to cluster j
                if Cj_ind == Ci_ind
                    change = 0
                else

                    Znew = copy(Z)
                    Znew[i] = Cj_ind
                    change = modularity(H,Znew,Ω) - modularity(H,Z,Ω)
                end

                # Check if this is currently the best possible greedy move to make
                if change > BestImprove
                    BestImprove = change
                    BestC = Cj_ind
                    improving = true
                end
            end

            # Move i to the best new cluster, only if it strictly improves modularity
            if BestImprove > 1e-8
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


function HyperLouvain(H::hypergraph,kmax::Int64,Ω,maxits::Int64=100,bigInt::Bool=true)
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
    println("")
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
            println("Louvain Iteration $iter")
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

                    # Znew = copy(Z)
                    # Znew[i] = Cj_ind
                    # voldiff = -second_term_eval(H, Znew, Ω; bigInt = bigInt)+second_term_eval(H, Z, Ω; bigInt = bigInt)

                    voldiff, ΔV, Δμ, ΔM = compute_voldiff(V, μ, M, i, Cj_ind, H.D, Z, C)
                    cdiff = CutDiff(Hyp,w,node2edges,Z,i,Cj_ind,Ω)
                    change =  cdiff + voldiff
                end

                # Check if this is currently the best possible greedy move to make
                if change > BestImprove
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


function compute_voldiff(V::Array, μ::Array, M::Dict,i::Int64, t::Int64, D::Vector{Int64}, Z::Vector{Int64},C::Dict)

    # increments due to proposal
    ΔV, Δμ, ΔM = increments(V, μ, M, i, t, D, Z)

    # new proposed quantities
    V_prop, μ_prop, M_prop = addIncrements(V, μ, M, ΔV, Δμ, ΔM)

    vol = 0
    for p in keys(M)
        vol += Ω(p,mode = "partition")*M[p]*C[p]
    end
    vol_prop = 0
    for p in keys(M_prop)
        vol_prop += Ω(p,mode = "partition")*M_prop[p]*C[p]
    end
    voldiff = vol-vol_prop

    return voldiff, ΔV, Δμ, ΔM
end
