using StatsBase
using Combinatorics

include("HSBM.jl")

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
    for i = 1:r, j = 1:i, p in partitions(i, j)
        N[p] = evalSumPV(p, Z, D)
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

function computeMoments(Z::Array, D::Array, r, ℓ=0)
    if ℓ==0
        ℓ=maximum(Z)
    end

    V = [sum(D.*(Z .== j)) for j in 1:ℓ]
    μ = [sum(V.^i) for i = 1:r]
    return(V, μ)
end

function evalConstants(r)
    C = Dict{Array{Integer, 1}, Integer}()
    for i = 1:r, j = 1:i, p in partitions(i, j)
        orderCorrection = prod([factorial(c) for c in values(countmap(p))])
        C[p] = multinomial(p...) ÷ orderCorrection
    end
    return(C)
end;

function evalSums(Z::Array, D::Array, r; constants=true, ℓ=0, bigInt=true)
    """
    Z: an Array of integer group labels
    D: an Array of degrees
    r: the largest hyperedge size to compute
    """

    if bigInt
        D = convert(Array{BigInt,1}, D)
        Z = convert(Array{BigInt,1}, Z)
    end

    V, μ = computeMoments(Z, D, r, ℓ)
    
    M = Dict{Array{Integer, 1}, Integer}()

    for i = 1:r, j = 1:i, p in partitions(i, j)
        M[p] = μ[p[end]]*get(M, p[1:(end-1)], 1) - correctOvercounting(M,p)
    end

    if constants
        C = evalConstants(r)
        for p in keys(M)
            M[p] *= C[p]
        end
    end

    return V, μ, M
end


# Extra methods for evalSums for working with degree dicts and the custom hypergraph class. 

function evalSums(Z::Array, D::Dict, r, ℓ=0, bigInt=true)
    d = [D[i] for i in 1:length(Z)]
    return evalSums(Z, d, r; ℓ=ℓ, bigInt=bigInt)
end

function evalSums(Z::Array, H::hypergraph, ℓ=0, bigInt=true)
    r = maximum(keys(H.E))
    return evalSums(Z, H.D, r, ℓ, bigInt)
end

function momentIncrements(V, μ, i, t, D, Z)
    """
    update V and μ in place by moving node i 
    from its current group to group t
    while returning the increment in μ for subsequent computation
    """
        
    ΔV = zero(V)
    # update V
    s      = Z[i];
    ΔV[t] += D[i];
    ΔV[s] -= D[i];
    
    r = length(μ)
    Δμ = [(V[s] + ΔV[s])^j - V[s]^j + (V[t] + ΔV[t])^j - V[t]^j for j in 1:r]

    return ΔV, Δμ
end

function increments(V, μ, M, i, t, D, Z)
    
    ΔV, Δμ = momentIncrements(V, μ, i, t, D, Z)

    # compute increments in M using recursion formula from notes
    ΔM = Dict{Array{Integer, 1}, Integer}()
    for i = 1:r, j = 1:i, p in partitions(i, j)
        ΔM[p] = Δμ[p[end]]*get(M, p[1:(end-1)], 1) + μ[p[end]]*get(ΔM, p[1:(end-1)], 0) + Δμ[p[end]]*get(ΔM, p[1:(end-1)], 0) - correctOvercounting(ΔM,p)
    end

    return(ΔV, Δμ, ΔM)
end

function addIncrements(V, μ, M, ΔV, Δμ, ΔM)
    M̃ = Dict(p => M[p] + ΔM[p] for p in keys(M))
    return(V + ΔV, μ + Δμ, M̃)
end

function second_term_eval(H::hypergraph, Z::Array{Int64, 1}, ℓ::Int64, Ω, bigInt=true)
    """
    Naive implementation, computes sums from scratch. 
    H: hypergraph
    Z: array storing cluster indices; c[i] is the cluster node i is in
    ℓ: maximum hyperedges size in H
    Ω: group interation function (e.g., planted partition). Needs to have a mode argument which, when set to value "partition", will cause evaluation on partition vectors rather than label vectors. 
    """

    obj = 0
    V, μ, M = evalSums(Z, H, ℓ, bigInt)
    for p in keys(M)
        obj += Ω(p, mode = "partition")*M[p]
    end
    return obj
end