"""
Functions for estimating an intensity function Ω from observed data. Throughout, H is an object of type hypergraph and Z is a label vector of the same length as H.D. The entries of Z are integers between 1 and kmax, where kmax is the number of clusters. Z is assumed to contain at least one entry with value k for each k ∈ [kmax]. 
"""

function estimateΩEmpirically(H, Z; pseudocount=0.0, aggregator=identity, bigNums = true)
    """
    Estimate an intensity function Ω directly from data, with no parameterization. 
    The intensity function is assumed to be symmetric under permutations of node labels. 
    Equivalently, it is a function of the partition vector p alone. 

    # Arguments
    - H::hypergraph, a hypergraph
    - Z::Vector{Int4}, a label vector
    - min_val::Float64, the pseudocount used when no hyperedges with a given specified feature key are observed. Altering the default is not recommended. 
    - aggregator::function, an aggregation function that maps a partition vector p to a feature vector summarizing p. 
    - bigNums::bool, whether to use bigInt and bigFloat in computation. Altering the default is not recommended except on exceptionally small data instances and with great caution. 

    # Returns

    - Ω::IntensityFunction, giving the maximum likelihood estimates of the intensity parameter for each possible value of the supplied aggregator.  
    """

    ℓ = maximum(keys(H.E))  # size of largest hyperedge
    S = evalSums(Z,H,maximum(Z),true)[3]
    
    C = Dict(p => 0.0 for p in partitionsUpTo(ℓ))
    for k in keys(H.E)
        Ek = H.E[k]
        for e in keys(Ek)
            p = partitionize(Z[e])
            C[p] = get(C, p, 0) + Ek[e]
        end
    end
    
    if bigNums
        C_agg = Dict(aggregator(p) => big(0.0) for p in keys(S))
        S_agg = Dict(aggregator(p) => big(0.0) for p in keys(S))
    else
        C_agg = Dict(aggregator(p) => 0.0 for p in keys(S))
        S_agg = Dict(aggregator(p) => 0.0 for p in keys(S))
    end
        
    for p in keys(S)
        agg_key = aggregator(p)
        if bigNums
            C_agg[agg_key] += get(C, p, big(pseudocount))
            S_agg[agg_key] += get(S, p, big(pseudocount))
        else
            C_agg[agg_key] += get(C, p, pseudocount)
            S_agg[agg_key] += get(S, p, pseudocount)
        end
    end
    
    ω̂ = Dict(a => abs(C_agg[a] / S_agg[a]) for a in keys(C_agg))
            
    return empiricalIntensityFunction((p, α) -> ω̂[p], ℓ, aggregator)
end

function learnParameters(H, Z, Ω, α0; max_iters = 10, verbose = false, tol = 1e-2)    
    """
    Estimate the parameter vector α for a parameterized IntensityFunction Ω from a hypergraph H and label vector Z. 
    Currently uses a very simple coordinatewise method, and is quite slow for larger data. Tuning amin and amax may be necessary to achieve appropriate results. 

    # Arguments

    - H::hypergraph, a hypergraph
    - Z::Vector{Int4}, a label vector
    - Ω::IntensityFunction, the intensity function parameterized by α
    - α0::Vector{Float64}, the initial guess for α. 
    - n_iters::Int64, the number of out iterations in coordinate ascent to perform
    - amin::Float64 and amax::Float64, box constraints on the entries of α

    # Returns

    α: Vector{Float64}, the parameter vector obtained via repeated optimization of the likelihood objective. 
    """
    kmax = length(α0) ÷ 2

    obj = formObjective(H, Z, Ω)
        
    function G(a, α, k)
        α_ = copy(α)
        α_[k] = a[1]
        return obj(α_)
    end

    α = copy(α0);
    
    for i = 1:max_iters
        old_val = obj(α)
        for k in vcat((kmax+1):2kmax, 1:kmax)
            res = Optim.optimize(a -> G(a, α, k), [α[k]], Optim.LBFGS()) # very slow and simple -- no gradient information
            α[k] = Optim.minimizer(res)[1]
        end
        new_val = obj(α)
        if abs((new_val - old_val)) < tol
            return α
        end
        if verbose
            println("Q = $Float64(modularity(H, Z, Ω;α = α))")
        end
    end
    return(α)
end


function formObjectives(H, Z, Ω)
    ℓ       = maximum(Z)
    C       = evalCuts(H,Z,Ω)
    V, μ, S = evalSums(Z,H,ℓ,true)
    Ŝ       = aggregateSums(S,Ω)
    kmax    = maximum(keys(H.E))

    function obj(α, k)
        val = 0.0
        for p in keys Ŝ
            if sum(p) == k
                Op   = Ω.ω(p, α)
                val += get(C, p, 0)*log(Op) - Ŝ[p]*Op
            end
        end
        return - val
    end
    return [α → obj(α, k) for k in 1:kmax]
end



function formObjective(H, Z, Ω)
    ℓ       = maximum(Z)
    C       = evalCuts(H,Z,Ω)
    V, μ, S = evalSums(Z,H,ℓ,true)
    Ŝ       = aggregateSums(S,Ω)
    
    function objective(α)
        obj = 0.0
        for p in keys(Ŝ)
            Op   = Ω.ω(p, α)
            obj += get(C, p, 0)*log(Op) - Ŝ[p]*Op
        end
        return -obj # sign is for minimization
    end
    return objective
end