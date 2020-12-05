
function modularity(H::hypergraph, Z::Array{<:Integer, 1}, Ω::IntensityFunction; α, bigInt::Bool=true)
    """
    Compute the modularity of a partition Z in a hypergraph H with interaction function Ω.
    H::hypergraph: the input hypergraph
    Z::array{Int64, 1}: array of group labels.
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended.
    return: Q::Float, the modularity term in the HSBM likelihood for H, Z, and Ω.
    """

    cut = first_term_eval(H, Z, Ω; α=α)
    vol = second_term_eval(H, Z, Ω; α=α, bigInt = bigInt)

    return cut - vol
end

function logLikelihood(H::hypergraph, Z::Array{<:Integer, 1}, Ω::IntensityFunction; α, bigInt::Bool=true)
    """
    Compute the HSBM log-likelihood of a partition Z in a hypergraph H with interaction function Ω.
    H::hypergraph: the input hypergraph
    Z::array{Int64, 1}: array of group labels.
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended.
    return: L::Float, the log-likelihood of H with parameters Z and Ω under the HSBM model.
    """
    Q = modularity(H, Z, Ω; α=α, bigInt=bigInt)

    D = computeDegrees(H)

    K, C = 0, 0

    logD = log.(D)

    kmax = maximum(keys(H.E))

    for ℓ = 1:kmax
        if haskey(H.E, ℓ)
            El = H.E[ℓ]
            for edge in keys(El)
                c = counting_coefficient(edge)
                if c < 0 print("woops") end
                weight = El[edge]
                K += weight*sum(logD[edge])
                C += weight*log(c) - log(factorial(weight)) # can maybe improve on performance here?
            end
        end
    end
    return Q, K, C
end
