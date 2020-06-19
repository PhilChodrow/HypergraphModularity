using Combinatorics
include("omega.jl")
include("HSBM.jl")
include("utils.jl")
include("hyper_format.jl")

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

function first_term_v2(H::Vector{Vector{Int64}},w::Array{Float64,1},c::Array{Int64,1}, Ω)
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

function CutDiff(H::Vector{Vector{Int64}},w::Array{Float64,1},node2edges::Vector{Vector{Int64}},c::Array{Int64,1}, I::Int64,J::Int64,Ω)
    """
    CutDiff: Compute change in the first term of the modularity function
    resulting from moving a node I to cluster J.

    This code is sloppy and I still think there's a faster way, but this is a step.
    """
    orig = c[I]
    obj1 = 0
    He = node2edges[I]
    for i = He
        edge = H[i]
        l = length(edge)
        clus_e = c[edge]    # set of clusters
        p = partitionize(clus_e)
        om_z = Ω(p; mode="partition", k = l)
        obj1 += w[i]*log(om_z)
    end

    obj2 = 0
    c[I] = J
    for i = He
        edge = H[i]
        l = length(edge)
        clus_e = c[edge]    # set of clusters
        p = partitionize(clus_e)
        om_z = Ω(p; mode="partition", k = l)
        obj2 += w[i]*log(om_z)
    end
    c[I] = orig
    return obj2-obj1
end


function NaiveCutDiff(H::hypergraph,c::Array{Int64,1},I::Int64, J::Int64, Ω)
    """
    NaiveCutDiff: Compute change in the first term of the modularity function
    resulting from moving a node I to cluster J.

    Naive edition: just call first_term_eval twice
    """
    orig = c[I]
    obj1 = first_term_eval(H,c,Ω)
    c[I] = J
    obj2 = first_term_eval(H,c,Ω)
    c[I] = orig
    return obj2 - obj1
end


function NaiveCutDiff2(HyperedgeList::Vector{Vector{Int64}},w::Array{Float64,1},node2edges::Vector{Vector{Int64}},c::Array{Int64,1},I::Int64,J::Int64, Ω)
    """
    NaiveCutDiff2: Compute change in the first term of the modularity function
    resulting from moving a node I to cluster J.

    Naive edition 2: just call first_term_v2 twice
    """
    orig = c[I]
    obj1 = first_term_v2(HyperedgeList,w,c,Ω)
    c[I] = J
    obj2 = first_term_v2(HyperedgeList,w,c,Ω)
    c[I] = orig
    return obj2 - obj1
end
