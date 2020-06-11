include("utils.jl")
include("cut.jl")
include("vol.jl")

function modularity(H::hypergraph, Z::Array{Int64, 1}, kmax::Int64, Ω;  bigInt=true)

    cut = first_term_eval(H, Z, Ω, kmax)
    vol = second_term_eval(H, Z, Ω, kmax, bigInt)
    return cut - vol
end

function L(H, Z, Ω, kmax, bigInt=true)
    Q = modularity(H, Z, kmax, Ω; bigInt=bigInt)   

    D = computeDegrees(H)

    K = 0
    C = 0

    logD = log.(D)

    for ℓ = 1:kmax
        El = H.E[ℓ]
        for edge in keys(El)
            weight = El[edge]
            K += weight*sum(logD[edge])
            C -= log(factorial(weight)) # can maybe improve on performance here?
        end
    end
    return Q, K, C
end

function logLikelihood(H::hypergraph, Z::Array{Int64,1}, Ω::Any) 
    """
    Given a hypergraph, return the HSBM likelihood using labels Z, degree parameters ϑ, and group intensities Ω.
    NOTE: this is a VERY slow function that should be spead up by orders of magnitude when Ω falls into important special cases

    BUGGED: needs combinatorial weighting to account for different numbers of passes over different combinations. 
    """
    n = length(Z)
    L = 0

    for k in keys(H.E)  
        T = with_replacement_combinations(1:n, k) 
        Ek = H.E[k]   
        for S in T

            z = Z[S]
            p = partitionize(z)
            factor = prod([factorial(c) for c in p])
            factor *= multinomial(p...)
            θ = H.D[S]
            X = Poisson(prod(θ)*Ω(z; mode="group"))
            
            m = get(Ek, S, 0)
            L += log(pdf(X, m))
        end
    end
    return(L)
end