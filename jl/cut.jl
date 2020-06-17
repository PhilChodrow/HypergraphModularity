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
    
            clus_e = c[edge]    # set of clusters
            weight = El[edge]
            obj   += weight*log(Ω(clus_e; mode="group"))
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
    e.g., given hyperedge size function fk and partition function fp, can use:
        ff = p->fp(p)*fk(sum(p))
        Om = build_omega(kmax,ff) # kmax = max hyperedge size
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


# """
# Rather than recomputing the penalty associated with each type of partition
# vector every time we see it, we can store penalties.

# Om[i][p] will return the Ω value for partition vector p for hyperedge size i.

# kmax = maximum hyperedge size
# fp = function for penalty associated with this partition
# """
# function build_omega(kmax,fp)

#     # each hyperedge size holds its own set of penalties
#     Om = Vector{Dict}()
#     for i = 1:kmax
#         ioms = Dict()

#         # generate all integer partitions of i
#         for part in partitions(i)
#             ioms[part] = fp(part)
#         end
#         push!(Om,ioms)
#     end
#     return Om
# end

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
