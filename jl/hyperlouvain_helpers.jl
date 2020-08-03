function partitionize(a::Array{Int64,1})
    """
    For a given integer vector a, return the partition corresponding to that
    vector. Useful for both counting corrections when sampling and computing
    likelihoods, and when computing partition-based values of Ω.

    This is the fastest version I could come up with.
    """
    k = length(a)
    d = Dict{Int64,Int64}()
    d[a[1]] = 1     # dictionary from cluster index to id in vector v
    next = 2
    v = [1]         # number of times the cluster index shows up
    for i = 2:k
        ind = get(d, a[i], next)
        if ind == next
            push!(v,1)
            next += 1
        else
            v[ind] += 1
        end
    end
    return sort(v,rev = true)
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


function CutDiff_OneNode(H::Vector{Vector{Int64}},w::Array{Float64,1},node2edges::Vector{Vector{Int64}},c::Array{Int64,1}, I::Int64,J::Int64,edge2penalty::Vector{Float64},Ω; α)
    """
    CutDiff: Compute change in the first term of the modularity function
    resulting from moving a node I to cluster J.
    """
    orig = c[I]
    He = node2edges[I]
    Δpen = Dict{Int64,Float64}()    # Change in penalty for moving a node

    obj = 0
    c[I] = J
    for i = He
        edge = H[i]
        clus_e = c[edge]    # set of clusters
        p = partitionize(clus_e)
        om_z = Ω(p; α=α, mode="partition")
        Δpen[i] = w[i]*log(om_z)   # new penalty if this move is made
        obj += w[i]*log(om_z)-edge2penalty[i]

    end
    c[I] = orig

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
    Hnew = sparse(J,I,ones(length(I)),m,sn)
    Neighbs = NeighborList(Hnew, sparse(Hnew'))
    return Neighbs
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
