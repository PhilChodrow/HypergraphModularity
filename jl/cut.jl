using Combinatorics
include("omega.jl")
include("HSBM.jl")
include("utils.jl")
include("hyper_format.jl")


function first_term_eval(H::hypergraph,c::Array{Int64,1}, Ω; α)
    """
    Not optimized, goal is to make this as quick and easy as
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
            obj += El[edge]*log(Ω(p; α=α, mode="partition"))
        end
    end
    return obj
end

function first_term_v2(H::Vector{Vector{Int64}},w::Array{Float64,1},c::Array{Int64,1}, Ω; α)
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
        clus_e = c[edge]    # set of clusters

        p = partitionize(clus_e)

        om_z = Ω(p; α=α, mode="partition")

        obj += w[i]*log(om_z)
    end
    return obj
end


function CutDiff(H::Vector{Vector{Int64}},w::Array{Float64,1},node2edges::Vector{Vector{Int64}},c::Array{Int64,1}, I::Int64,J::Int64,Ω; α)
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
        clus_e = c[edge]    # set of clusters
        p = partitionize(clus_e)
        om_z = Ω(p;α=α, mode="partition")
        obj1 += w[i]*log(om_z)
    end

    obj2 = 0
    c[I] = J
    for i = He
        edge = H[i]
        clus_e = c[edge]    # set of clusters
        p = partitionize(clus_e)
        om_z = Ω(p; α=α, mode="partition")
        obj2 += w[i]*log(om_z)
    end
    c[I] = orig
    return obj2-obj1
end


function NaiveCutDiff(H::hypergraph,c::Array{Int64,1},I::Int64, J::Int64, Ω; α)
    """
    NaiveCutDiff: Compute change in the first term of the modularity function
    resulting from moving a node I to cluster J.

    Naive edition: just call first_term_eval twice
    """
    orig = c[I]
    obj1 = first_term_eval(H,c,Ω; α=α)
    c[I] = J
    obj2 = first_term_eval(H,c,Ω; α=α)
    c[I] = orig
    return obj2 - obj1
end


function NaiveCutDiff2(HyperedgeList::Vector{Vector{Int64}},w::Array{Float64,1},node2edges::Vector{Vector{Int64}},c::Array{Int64,1},I::Int64,J::Int64, Ω; α)
    """
    NaiveCutDiff2: Compute change in the first term of the modularity function
    resulting from moving a node I to cluster J.

    Naive edition 2: just call first_term_v2 twice
    """
    orig = c[I]
    obj1 = first_term_v2(HyperedgeList,w,c,Ω; α=α)
    c[I] = J
    obj2 = first_term_v2(HyperedgeList,w,c,Ω; α=α)
    c[I] = orig
    return obj2 - obj1
end

function evalCuts(Z::Array{Int64,1}, H::hypergraph)
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

function first_term_v3(Z::Array{Int64,1},H::hypergraph, Ω; α)
    """
    This function may be slightly more efficient than v2 due to the fact that we multiply by log Ω fewer times. It is here primarily as a check on evalCuts, which has independent significance in the context of learning parameterized versions of Ω. 
    """
    C = evalCuts(Z, H)
    sum(C[p]*log(Ω(p; α=α, mode="partition")) for p in keys(C))
end