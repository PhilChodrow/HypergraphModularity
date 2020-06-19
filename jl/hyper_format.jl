# These are some useful functions for converting back and forth between
# different ways to store hypergraphs.

using SparseArrays
using LinearAlgebra

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
    Hyperedges = Vector{Vector{Int64}}()
    n,m = size(H)

    for i = 1:m
        startedge = ci[i]
        endedge = ci[i+1]-1
        edge = rp[startedge:endedge]
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
    NV: Also yes, historically I have found arrays to be more forgiving than dictionaries in Julia, but I also wouldn't be surprised if I find a counterexample to that at some point. In theory the dictionary should eventually be better.
    """
    Hyperedges = Vector{Vector{Int64}}()
    weights = Vector{Int64}()
    for key in keys(H.E)
        ed = collect(keys(H.E[key]))
        wt = collect(values(H.E[key]))
        append!(weights,wt)
        append!(Hyperedges,ed)
    end
    return Hyperedges, weights
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
