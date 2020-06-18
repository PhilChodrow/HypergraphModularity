using Combinatorics
include("omega.jl")
include("HSBM.jl")
include("utils.jl")

function first_term_eval(H::hypergraph,c::Array{Int64,1}, Ω)

    """
    First version: not optimized, goal is to make this as quick and easy as
    possible using existing code.
    H: hypergraph
    c: array storing cluster indices; c[i] is the cluster node i is in
    kmax: maximum hyperedges size in H
    Ω: group interation function (e.g., planted partition)
    """

    kmin, kmax = minimum(keys(H.E)), maximum(keys(H.E))

    obj = 0
    for l = kmin:kmax
        El = H.E[l]
        for edge in keys(El)
            p = partitionize(c[edge])
            obj += El[edge]*log(Ω(p; mode="partition", k = l))
        end
    end
    return obj
end

function first_term_v2(H::Vector{Vector{Int64}},w::Array{Int64,1},c::Array{Int64,1}, Ω)
    """
    Second version: store penalties first,
    and more importantly, a faster way to compute Ω(z_e)
    H: hypergraph just stored as an edge list
    c: array storing cluster indices; c[i] is the cluster node i is in
    kmax: maximum hyperedges size in H
    Ω: Interaction function. Must have been built with option by_size = true using buildΩ().  
    """
    obj = 0

    for i = 1:length(w)
        
        edge = H[i]
        l = length(edge)
        clus_e = c[edge]    # set of clusters

        p = partitionize(clus_e)

        om_z = Ω(p; mode="partition", k = l) 

        obj += w[i]*log(om_z)
    end
    return obj
end

"""
Convert a hypergraph from old to new format.
PC: I find that first_term_v2 is around twice as fast as first_term_eval -- if that turns out to make enough of a difference, then maybe this should actually be our core data structure? 
PC: If I understand correctly, this is advantageous to do because indexing Arrays is faster than Dict lookups?
"""
function hyperedge_formatting(H::hypergraph)

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