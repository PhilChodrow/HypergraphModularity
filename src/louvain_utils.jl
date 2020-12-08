"""
The way the Louvain code works, you start with n clusters (singletons), and as nodes
move, many of these become empty. At certain points in the algorithm, it is
helpful to renumber the cluster IDs, so that they go from 1 to M where
M = number of non-empty clusters.

e.g. c = [2 3 5 3 9] ---> c = [1 2 3 2 4]

This version does this with both the Clusters array of arrays AND
the cluster indicator vector c:

c[i] = (integer) cluster ID that node i belongs to
Clusters[j] = (integer array) node IDs for nodes in cluster j
"""
function renumber(c::Vector{Int64},Clusters::Vector{Vector{Int}})

    n = length(c)
    map = sort(unique(c))     # map from new cluster ID to old cluster ID
    old2new = Dict()    # map from old cluster ID to new cluster ID
    for i = 1:length(map)
        old2new[map[i]] = i
    end
    cnew = zeros(Int64,n)

    Clusters = Clusters[map]

    # Rename the clusters
    for i = 1:n
        # newClus = findfirst(x->x == c[i],map)
        newClus = old2new[c[i]]
        cnew[i] = newClus
        push!(Clusters[newClus],i)
    end

    return cnew, Clusters

end


"""
See above function, this function doesn't care about the Clusters array
"""
function renumber(c::Vector{Int64})

    n = length(c)
    map = sort(unique(c))     # map from new cluster ID to old cluster ID
    old2new = Dict()          # map from old cluster ID to new cluster ID
    for i = 1:length(map)
        old2new[map[i]] = i
    end
    cnew = zeros(Int64,n)

    # Rename the clusters
    for i = 1:n
        cnew[i] = old2new[c[i]]
    end

    return cnew

end

# These are some useful functions for converting back and forth between
# different ways to store hypergraphs.

function NeighborList(He2n::SparseArrays.SparseMatrixCSC{Float64,Int64},Hn2e::SparseArrays.SparseMatrixCSC{Float64,Int64})
    """
    NeighborList: given hypergraph (edge x node) indicence and (node x edge)
    incidence, return a list of neighbors for each node.
    """
    m,n = size(He2n)
    rp = He2n.rowval
    ci = He2n.colptr
    Neighbs = Vector{Vector{Int64}}()
    for i = 1:n
        # Get list of edges that node i is in
        Edges = rp[ci[i]:ci[i+1]-1]

        # build a set of neighboring nodes
        NeighbSet = Vector{Int64}()
        for j = 1:length(Edges)
            ej = Edges[j]               # For each node ej adjacent to I...
            nodes = getnodes(Hn2e,ej)   # ...get nodes in edge ej
            append!(NeighbSet,nodes)
        end
        push!(Neighbs, sort(setdiff(unique(NeighbSet),i)))
    end

    return Neighbs
end


function NeighborList(node2edge::Vector{Vector{Int64}},edge2node::Vector{Vector{Int64}})
    """
    NeighborList: given edge to node list and node to edge list, compute
    neighbors
    """
    n = length(node2edge)
    m = length(edge2node)
    Neighbs = Vector{Vector{Int64}}()
    for i = 1:n
        # Get list of edges that node i is in
        Edges = node2edge[i]

        # build a set of neighboring nodes
        NeighbSet = Vector{Int64}()
        for j = 1:length(Edges)
            ej = Edges[j]               # For each node ej adjacent to I...
            nodes = edge2node[ej]       # ...get nodes in edge ej
            append!(NeighbSet,nodes)
        end
        push!(Neighbs, sort(setdiff(unique(NeighbSet),i)))
    end

    return Neighbs
end

function getedges(He2n::SparseArrays.SparseMatrixCSC{Float64,Int64},I::Int64)
    """
    Given the edge-by-node incidence matrix He2n, return the set of hyperedges
    that node I is in.
    """
    first = He2n.colptr[I]
    last = He2n.colptr[I+1]-1
    edges = He2n.rowval[first:last]

    return edges
end

function getnodes(Hn2e::SparseArrays.SparseMatrixCSC{Float64,Int64},J::Int64)
    """
    Given the node-by-edge incidence matrix Hn2e, return the set of nodes
    contained in hyperedge J
    """
    first = Hn2e.colptr[J]
    last = Hn2e.colptr[J+1]-1
    nodes = Hn2e.rowval[first:last]
    mult = Hn2e.nzval[first:last]
    edge = Vector{Int64}()
    # need to adjust for multiplicities
    for t = 1:length(nodes)
        node = nodes[t]
        for k = 1:mult[t]
            push!(edge,node)
        end
    end

    return edge
end

function incidence2elist(H::SparseArrays.SparseMatrixCSC{Float64,Int64},nodelist::Bool=false)
    """
    This converts a hypergraph in incidence matrix form to hyperedge list form.
    Incidence matrix form:  H[e,u] = 1  iff node u is in hyperedge e
    Edgelist form: Hyperedges[j] = array of nodes in hyperedge j

    nodelist == true: means you actually want to get a map from node IDs to hyperedge ids
    """
    if ~nodelist
        # unless you want the node2edge map, transpose first
        H = SparseArrays.sparse(H')
    end
    rp = H.rowval
    ci = H.colptr
    nz = H.nzval
    Hyperedges = Vector{Vector{Int64}}()
    n,m = size(H)

    for i = 1:m
        startedge = ci[i]
        endedge = ci[i+1]-1
        nodes = rp[startedge:endedge]
        mult = nz[startedge:endedge]
        edge = Vector{Int64}()
        # need to adjust for multiplicities
        for t = 1:length(nodes)
            node = nodes[t]
            for k = 1:mult[t]
                push!(edge,node)
            end
        end

        push!(Hyperedges,edge)
    end
    return Hyperedges
end

function elist2incidence(Hyperedges::Vector{Vector{Int64}}, N::Int64)
    """
    Converts a hyperedge list into a binary incidence matrix for the hypergraph.

    This is the exact inverse of incidence2elist
    """
    U = Vector{Int64}()
    E = Vector{Int64}()
    M = length(Hyperedges)
    for enum = 1:length(Hyperedges)
        e = Hyperedges[enum]
        for node in e
            push!(U,node)
            push!(E,enum)
        end
    end

    H = SparseArrays.sparse(E,U,ones(length(U)),M,N)
    return H
end

function hyperedge_formatting(H::hypergraph)
    """
    Convert a hypergraph from dictionary-based format (type 'hypergraph')
    to Hyperedge list with weights vector
    """
    Hyperedges = Vector{Vector{Int64}}()
    weights = Vector{Float64}()
    for key in keys(H.E)
        ed = collect(keys(H.E[key]))
        wt = collect(values(H.E[key]))
        append!(weights,wt)
        append!(Hyperedges,ed)
    end
    return Hyperedges, weights
end


function hypergraph2incidence(H::hypergraph)
    """
    Convert dictionary hypergraph format to binary edge-by-node incidence matrix
    """
    Hyperedges = Vector{Vector{Int64}}()
    weights = Vector{Float64}()
    for key in keys(H.E)
        ed = collect(keys(H.E[key]))
        wt = collect(values(H.E[key]))
        append!(weights,wt)
        append!(Hyperedges,ed)
    end
    N = length(H.D)
    He2n = elist2incidence(Hyperedges,N)
    return He2n, weights
end

function EdgeMap(H::hypergraph)
    """
    Given a hypergraph, extract a map from node indices to hyperedge indices
    """
    n = length(H.D)
    Hyp, weights = hyperedge_formatting(H)
    Hin = elist2incidence(Hyp,n)

    return incidence2elist(Hin,true)
end
