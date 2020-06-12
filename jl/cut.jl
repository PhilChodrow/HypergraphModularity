using Combinatorics
include("omega.jl")
include("HSBM.jl")
include("utils.jl")


function first_term_eval(H::hypergraph,c::Array{Int64,1},kmax::Int64,kmin::Int64, Ω)

    """
    First version: not optimized, goal is to make this as quick and easy as
    possible using existing code.
    H: hypergraph
    c: array storing cluster indices; c[i] is the cluster node i is in
    kmax: maximum hyperedges size in H
    Ω: group interation function (e.g., planted partition)
    """
    obj = 0
    # Is there any reason to start with l = 1 sized hyperedges?
    for l = kmin:kmax
        El = H.E[l]
        lfac = factorial(l)
        for edge in keys(El)
            perms = count_coefficient(edge)
            clus_e = c[edge]    # set of clusters
            weight = El[edge]
            obj   += perms*weight*log(Ω(clus_e; mode="group"))
        end
    end
    return obj
end

function first_term_v2(H::Vector{Vector{Int64}},w::Array{Float64,1},c::Array{Int64,1},kmax::Int64, kmin::Int64,
    Om::Vector{Dict})
    """
    Second version: store penalties first,
    and more importantly, a faster way to compute Ω(z_e)
    H: hypergraph just stored as an edge list
    c: array storing cluster indices; c[i] is the cluster node i is in
    kmax: maximum hyperedges size in H
    Om: weights for group interation function
        Om[i][p] = Ω value for partition vector p, hyperedg size i.

    e.g., given hyperedge size function fk and partition function fp, can use:
        ff = p->fp(p)*fk(sum(p))
        Om = build_omega(kmax,ff) # kmax = max hyperedge size
    """
    obj = 0

    for i = 1:length(w)

        edge = H[i]
        l = length(edge)
        clus_e = c[edge]    # set of clusters

        p = cvec_2_pvec(clus_e,l)

        perms = counting_coefficient(edge)

        om_z = Om[l][p]

        obj += perms*w[i]*log(om_z)
    end
    return obj
end


"""
Rather than recomputing the penalty associated with each type of partition
vector every time we see it, we can store penalties.

Om[i][p] will return the Ω value for partition vector p for hyperedge size i.

kmax = maximum hyperedge size
fp = function for penalty associated with this partition
"""
function build_omega(kmax,fp)

    # each hyperedge size holds its own set of penalties
    Om = Vector{Dict}()
    for i = 1:kmax
        ioms = Dict()

        # generate all integer partitions of i
        for part in partitions(i)
            ioms[part] = fp(part)
        end
        push!(Om,ioms)
    end

    return Om
end

"""
cvec2pvec: cluster vector to partition vector.
Given an array of integers z, compute the integer partition vector p.
e.g. input z = [1 2 3 2 2] returns output p = [3 1 1].
"""
function cvec_2_pvec(a::Array{Int64,1},k::Int64)
    # we might be able to do this faster. Come back later.
    # p = countmap(a)
    # p = sort(collect(values(p)),rev=true)
    u = unique(a)
    d = Dict()
    lu = length(u)
    for j = 1:lu
        d[u[j]] = j
    end
    cnts = zeros(Int, lu)
    for i = 1:k
        @inbounds x = a[i]
        @inbounds cnts[d[x]] += 1
    end
    return sort(cnts,rev = true)
end


"""
Convert a hypergraph from old to new format.
"""
function hyperedge_formatting(H::hypergraph)

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
