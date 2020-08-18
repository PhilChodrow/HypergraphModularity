
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

# function parameterEstimateObjective(H::hypergraph, Z::Array{<:Integer, 1}, Ω; ℓ::Int64 = 0, bigInt::Bool=true)
#     """
#     An efficient way to compute the modularity objective for varying intensity function parameter α, useful for learning α from partitions.
#     Probably there is a better way to implement this via currying.
#     """

#     if ℓ == 0
#         ℓ = maximum(Z)
#     end

#     C = evalCuts(Z, H)
#     V, μ, M = evalSums(Z, H, ℓ, bigInt)

#     function objective(α)
#         obj = 0
#         for p in keys(M)
#             Op = Ω(p;α=α, mode="partition")
#             obj -= M[p]*Op
#             if p in keys(C)
#                 obj += C[p]*log(Op)
#             end
#         end
#         return -obj # for minimization
#     end
#     return objective
# end

function formObjective(H, Z, Ω)
    ℓ = maximum(Z)
    C       = evalCuts(Z,H)
    V, μ, S = evalSums(Z,H,ℓ,true);
    function objective(α)
        obj = 0.0
        for p in keys(S)
            Op   = Ω(p; α=α, mode="partition")
            obj += get(C, p, 0)*log(Op) - S[p]*Op
        end
        return -Float64(obj, RoundDown) # sign is for minimization
    end
    return objective
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
