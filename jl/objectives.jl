include("utils.jl")
include("cut.jl")
include("vol.jl")

function modularity(H::hypergraph, Z::Array{Int64, 1}, Ω; bigInt=true)

    kmax = maximum(keys(H.E))

    cut = first_term_eval(H, Z, Ω)
    vol = second_term_eval(H, Z, Ω; bigInt = bigInt)
    
    return cut - vol
end

function L(H, Z, Ω; bigInt=true)
    Q = modularity(H, Z, Ω; bigInt=bigInt)   

    D = computeDegrees(H)

    K, C = 0, 0
    
    logD = log.(D)

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
    Given a hypergraph, return the HSBM likelihood using labels Z, degree parameters ϑ, and group intensities Ω.
    NOTE: this is a VERY slow function that should be spead up by orders of magnitude when Ω falls into important special cases
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

# function logLikelihoodNaive(H::hypergraph, Z::Array{Int64,1}, Ω::Any) 
#     """
#     Given a hypergraph, return the HSBM likelihood using labels Z, degree parameters ϑ, and group intensities Ω.
#     This function just iterates over all possible tuples in all permutations -- used only for testing purposes. 
#     """
#     n = length(Z)
#     L = 0.0
#     D = 1.0*H.D

#     for k in keys(H.E)  
#         T = Iterators.product((1:n for i = 1:k)...)
#         Ek = H.E[k]  

#         for s in T
#             S = collect(s)
#             # S = sort(S)
#             c = counting_coefficient(S)

#             z = Z[S]
#             θ = D[S]

#             m = get(Ek, S, 0)
#             L += log(poisson_pdf(m, prod(θ)*Ω(z; mode="group")))
#         end
#     end
#     return(L)
# end