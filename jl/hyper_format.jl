# These are some useful functions for converting back and forth between
# different ways to store hypergraphs.

using SparseArrays
using LinearAlgebra

function NeighborList(He2n::SparseMatrixCSC{Float64,Int64},Hn2e::SparseMatrixCSC{Float64,Int64})
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

function getedges(He2n::SparseMatrixCSC{Float64,Int64},I::Int64)
    """
    Given the edge-by-node incidence matrix He2n, return the set of hyperedges
    that node I is in.
    """
    first = He2n.colptr[I]
    last = He2n.colptr[I+1]-1
    edges = He2n.rowval[first:last]

    return edges
end

function getnodes(Hn2e::SparseMatrixCSC{Float64,Int64},J::Int64)
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

function incidence2elist(H::SparseMatrixCSC{Float64,Int64},nodelist::Bool=false)
    """
    This converts a hypergraph in incidence matrix form to hyperedge list form.
    Incidence matrix form:  H[e,u] = 1  iff node u is in hyperedge e
    Edgelist form: Hyperedges[j] = array of nodes in hyperedge j

    nodelist == true: means you actually want to get a map from node IDs to hyperedge ids
    """
    if ~nodelist
        # unless you want the node2edge map, transpose first
        H = sparse(H')
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

    H = sparse(E,U,ones(length(U)),M,N)
    return H
end

function hyperedge_formatting(H::hypergraph)
    """
    Convert a hypergraph from old to new format.
    PC: I find that first_term_v2 is around twice as fast as first_term_eval -- if that turns out to make enough of a difference, then maybe this should actually be our core data structure?
    PC: If I understand correctly, this is advantageous to do because indexing Arrays is faster than Dict lookups?
    NV: Yes, might be worth considering making this the standard format, but probably not first priority at this point though.
    NV: Also yes, historically I have found arrays to be more forgiving than dictionaries in Julia
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
    Convert dictionary hypergraph format to binary edge x node incidence matrix
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
