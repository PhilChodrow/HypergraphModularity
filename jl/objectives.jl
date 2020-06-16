include("utils.jl")
include("cut.jl")
include("vol.jl")

function modularity(H::hypergraph, Z::Array{Int64, 1}, Ω; bigInt::Bool=true)
    """
    Compute the modularity of a partition Z in a hypergraph H with interaction function Ω. 
    H::hypergraph: the input hypergraph 
    Z::array{Int64, 1}: array of group labels.  
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended. 
    return: Q::Float, the modularity term in the HSBM likelihood for H, Z, and Ω. 
    """

    kmax = maximum(keys(H.E))

    cut = first_term_eval(H, Z, Ω)
    vol = second_term_eval(H, Z, Ω; bigInt = bigInt)
    
    return cut - vol
end

function L(H::hypergraph, Z::Array{Int64, 1}, Ω; bigInt::Bool=true)
    """
    Compute the HSBM log-likelihood of a partition Z in a hypergraph H with interaction function Ω. 
    H::hypergraph: the input hypergraph 
    Z::array{Int64, 1}: array of group labels.  
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended. 
    return: L::Float, the log-likelihood of H with parameters Z and Ω under the HSBM model. 
    """
    Q = modularity(H, Z, Ω; bigInt=bigInt)  
    

    # the remainder of the code focuses on computing the other terms in the likelihood. 


    # initialization prior to main loop
    D = computeDegrees(H)
    K, C = 0, 0
    
    logD = log.(D)

    # because of some convenient math simplifications, we are able to loop only through edges that exist, not over possible edges. 
    for ℓ = 1:kmax
        El = H.E[ℓ]
        for edge in keys(El)
            c = counting_coefficient(edge)
            weight = El[edge]
            K += weight*sum(logD[edge])
            C += weight*log(c) - log(factorial(weight)) # can maybe improve on performance here?
        end
    end
    return Q, K, C
end

function logLikelihood(H::hypergraph, Z::Array{Int64,1}, Ω::Any) 
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
    D = 1.0*H.D

    for k in keys(H.E)  
        T = with_replacement_combinations(1:n, k) 
        Ek = H.E[k]   
        for S in T

            z = Z[S]
            c = counting_coefficient(S)
            θ = D[S]
            
            m = get(Ek, S, 0)

            L += log(poisson_pdf(m, c*prod(θ)*Ω(z; mode="group")))

            K += m*(log(c) + sum(log.(θ)))
            V += c*prod(θ)*Ω(z; mode="group")
            R += log(factorial(m)) # think about this
            C += m*log(Ω(z; mode="group"))        
        end
    end
    return(L)
end