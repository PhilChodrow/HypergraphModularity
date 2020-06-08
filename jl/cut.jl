"""
Evaluating the first term of the modularity objective. This is
equivalent to a constant minus a hypergraph cut function.
"""

function first_term_eval(H::hypergraph,c::Array{Int64,1},kmax::Int64, Ω)
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
    for l = 1:kmax
        El = H.E[l]
        lfac = factorial(l)
        for edge in keys(El)
            clus_e = c[edge]    # set of clusters
            weight = El[edge]
            obj += lfac*weight*log(Ω(clus_e))
        end
    end
    return obj
end
