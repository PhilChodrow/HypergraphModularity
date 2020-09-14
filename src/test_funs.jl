"""
The functions in this script are not intended for use in computation.
Rather, they serve as naive and ground-truth computations of key quantities, which are used for unit testing in unit_tests.jl
"""

using StatsBase
using Combinatorics

# ------------------------------------------------------------------------------
# NAIVE VOLUME SUMS
# ------------------------------------------------------------------------------

# naive evaluation of the "volume" term (second term) in the modularity
# functional for partition-based interaction functions.
# It is not recommended to call these functions except on very small
# instances for testing purposes.

function evalSumNaive(p, Z, D)
    n = length(Z)
    r = sum(p)

    S = 0

    T = Iterators.product((1:n for i = 1:r)...)

    for R in T
        if partitionize(Z[collect(R)]) == p
            S += prod(D[collect(R)])
        end
    end
    return(S)
end

function evalSumNaive2(p, Z, D)
    n = length(Z)
    r = sum(p)
    s = 0

    T = with_replacement_combinations(1:n, r)

    for S in T
        if partitionize(Z[S]) == p
            s += prod(D[S])*counting_coefficient(S)
        end
    end
    return(s)
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

function test_sums(S, Z, D)
    n = length(Z)
    r = length(S)
    s1 = 0

    T = Iterators.product((1:n for i = 1:r)...)
    for q in T
        if sort(collect(q)) == S
            s1 += prod(D[S])
        end
    end

    s2 = 0
    T = with_replacement_combinations(1:n, r)
    for q in T
        if q == S
            s2 += prod(D[S]) * counting_coefficient(S)
        end
    end
    return(s1, s2)
end

# Sum-product-of-volumes representation.
# Much faster than the naive summation algorithms, but still very slow even on modest instances.
# It is again not recommended to call these functions directly except for testing purposes.

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


function naiveSecondTerm(H::hypergraph, Z::Vector{<:Integer}, Ω::IntensityFunction; α)
    """
    A naive evaluation of the second term of the modularity functional, for testing purposes.
    Looks like it is VERY bugged in some way that is strangely difficult to pin down.
    """

    n = length(H.D)
    obj = 0.0

    for k in keys(H.E)
        T =  Combinatorics.with_replacement_combinations(1:n, k)
        Ek = H.E[k]
        for S in T
            z = Z[S]
            c = counting_coefficient(S)
            d = H.D[S]
            p = Ω.P(z)

            obj += Ω.ω(p,α)*prod(d)*c
        end
    end
    return obj
end

function logLikelihoodNaive(H::hypergraph, Z::Array{<:Integer,1}, Ω::IntensityFunction, ϑ::Array{Float64,1} = zeros(1); α)
    """
    Compute the HSBM log-likelihood of a partition Z in a hypergraph H with interaction function Ω.
    This function is VERY slow and should generally only be used for testing purposes.
    H::hypergraph: the input hypergraph
    Z::array{Int64, 1}: array of group labels.
    Ω: group interaction function, as constructed by ΩFromDict(D)
    return: L::Float, the log-likelihood of H with parameters Z and Ω under the HSBM model.
    """
    n = length(Z)
    L = 0.0

    if ϑ == zeros(1)
        ϑ = 1.0*H.D
    end

    for k in keys(H.E)
        T = Combinatorics.with_replacement_combinations(1:n, k)
        Ek = H.E[k]
        for S in T

            z = Z[S]
            c = counting_coefficient(S)
            θ = ϑ[S]

            m = get(Ek, S, 0)

            p = Ω.P(z)
            L += log(poisson_pdf(m, c*prod(θ)*Ω.ω(p, α)))
        end
    end

    return(L)
end

function modularityNaive(H, Z, Ω; α)
    n = length(Z)
    Q = 0.0
    D = H.D

    for k in keys(H.E)
        T = Combinatorics.with_replacement_combinations(1:n, k)
        Ek = H.E[k]
        for S in T
            z = Z[S]
            c = counting_coefficient(S)
            θ = D[S]

            m = get(Ek, S, 0)
            p = Ω.P(z)
            om = Ω.ω(p, α)

            Q += m*log(om) - c*prod(θ)*om
        end
    end
    return Q
end
