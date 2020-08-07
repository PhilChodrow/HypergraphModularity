
function modularity(H::hypergraph, Z::Array{<:Integer, 1}, Ω; α, bigInt::Bool=true)
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

function parameterEstimateObjective(H::hypergraph, Z::Array{<:Integer, 1}, Ω; ℓ::Int64 = 0, bigInt::Bool=true)
    """
    An efficient way to compute the modularity objective for varying intensity function parameter α, useful for learning α from partitions.
    Probably there is a better way to implement this via currying.
    """

    if ℓ == 0
        ℓ = maximum(Z)
    end

    C = evalCuts(Z, H)
    V, μ, M = evalSums(Z, H, ℓ, bigInt)

    function objective(α)
        obj = 0
        for p in keys(M)
            Op = Ω(p;α=α, mode="partition")
            obj -= M[p]*Op
            if p in keys(C)
                obj += C[p]*log(Op)
            end
        end
        return -obj # for minimization
    end
    return objective
end



function L(H::hypergraph, Z::Array{<:Integer, 1}, Ω; α, bigInt::Bool=true)
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

function logLikelihood(H::hypergraph, Z::Array{<:Integer,1}, Ω::Any, ϑ::Array{Float64,1} = zeros(1); α)
    """
    Compute the HSBM log-likelihood of a partition Z in a hypergraph H with interaction function Ω.
    This function is VERY slow and should generally only be used for testing purposes.
    H::hypergraph: the input hypergraph
    Z::array{Int64, 1}: array of group labels.
    Ω: group interaction function, as constructed by ΩFromDict(D)
    return: L::Float, the log-likelihood of H with parameters Z and Ω under the HSBM model.
    """
    n = length(Z)
    L, C, V, K, R = 0.0, 0.0, 0.0, 0.0, 0.0

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

            L += log(poisson_pdf(m, c*prod(θ)*Ω(z; α = α, mode="group")))
        end
    end

    return(L)
end
