using StatsBase
using Combinatorics

"""
At the moment, the functions for evaluating sums in this module are strictly for the group-size partition intensity function and its relatives. 
"""

# ------------------------------------------------------------------------------
# NAIVE SUMMATION
# ------------------------------------------------------------------------------

# It is not recommended to call these functions except on very small instances for testing purposes. 

function evalSumNaive(p, Z, D)
    n = length(Z)
    r = sum(p)

    S = 0

    T = Iterators.product((1:n for i = 1:r)...)

    for R in T
        a = countmap(vec(Z[collect(R)]))
        a = -sort(-collect(values(a)))
        if a == p
            S += prod(D[collect(R)])
        end
    end
    return(S)
end

function evalSumsNaive(Z, D, r)
    N = Dict()
    for i = 1:r
        for j = 1:i # number of nonzero entries
            for p in partitions(i, j)
                N[p] = evalSumNaive(p, Z, D)
            end
        end
    end
    return(N)
end

# ------------------------------------------------------------------------------
# SUM-PRODUCT OF VOLUMES REPRESENTATION
# ------------------------------------------------------------------------------

# Much faster than the naive summation algorithms. 
# However, it is again not recommended to call these functions directly except for testing purposes. 

function evalSumPV(p, Z, D)
    n = length(Z)
    ℓ = sum(p)
    S = 0

    V = [sum([(Z[i] == s)*D[i] for i = 1:n]) for s = 1:maximum(Z)] # vector of volumes

    P = Iterators.product((1:maximum(Z) for i = 1:ℓ)...)
    for p_ in P
        a = countmap(vec(collect(p_)))
        a = -sort(-collect(values(a)))
        if a == p
          S += prod(V[collect(p_)])
        end
    end
    return(S)
end

function evalSumsPV(Z, D, r)
    N = Dict()
    for i = 1:r
        for j = 1:i # number of nonzero entries
            for p in partitions(i, j)
                N[p] = evalSumPV(p, Z, D)
            end
        end
    end
    return(N)
end



# ------------------------------------------------------------------------------
# FASTER RECURSIVE SUMMATION
# ------------------------------------------------------------------------------

# Compute a complete set of sums using a fast recursive algorithm based on a simple combinatorial relationship between the required sums. These are the only functions currently defined in this module which should be used in instances of larger than 10 nodes or so. 

function correctOvercounting(M::Dict, p::Array)
    """
    Utility function: second term in the recurrence in the notes
    """
    pk = p[end]
    S = 0
    for i = 1:length(p)-1
        p_ = copy(p)[1:(end-1)]
        p_[i] += pk
        S += M[-sort(-p_)]
    end
    return(S)
end

function computeμ(Z::Array, D::Array, r)
    V = [sum([D[i]*(Z[i] == j) for i in 1:length(Z)]) for j in 1:maximum(Z)]
    μ = [sum(V.^i) for i = 1:r]
    return(μ)
end

function evalSums(Z::Array, D::Array, r)
    """
    Z: an Array of integer group labels
    D: an Array of degrees
    r: the largest hyperedge size to compute
    """
    μ = computeμ(Z, D, r)
    
    M = Dict{Array{Integer, 1}, Integer}()

    for i = 1:r, j = 1:i
        for p in partitions(i, j)
            M[p] = μ[p[end]]*get(M, p[1:(end-1)], 1) - correctOvercounting(M,p)
        end
    end
    
    N = Dict{Array{Integer, 1}, Integer}()
    for p in keys(M)
        orderCorrection = 1
        counts = values(countmap(p))
        for c in counts
            orderCorrection *= factorial(c) 
        end
        N[p] = M[p] * multinomial(p...) ÷ orderCorrection
    end
    return(N)
end;

# Extra methods for evalSums for working with degree dicts and the custom hypergraph class. 

function evalSums(Z::Array, D::Dict, r)
    d = [D[i] for i in 1:length(Z)]
    return evalSums(Z, d, r)
end

function evalSums(Z::Array, H::hypergraph)
    r = maximum(keys(H.E))
    return evalSums(Z, H.D, r)
end