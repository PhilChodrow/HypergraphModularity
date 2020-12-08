function estimateΩEmpirically(H, Z; min_val=0.0, aggregator=identity, bigNums = true)
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
            C_agg[agg_key] += get(C, p, big(min_val))
            S_agg[agg_key] += get(S, p, big(min_val))
        else
            C_agg[agg_key] += get(C, p, min_val)
            S_agg[agg_key] += get(S, p, min_val)
        end
    end
    
    ω̂ = Dict(a => abs(C_agg[a] / S_agg[a]) for a in keys(C_agg))
            
    return empiricalIntensityFunction((p, α) -> ω̂[p], ℓ, aggregator)
end

function learnParameters(H, Z, Ω, α0; n_iters = 10, amin = 0, amax = 10)    
    
    kmax = length(α0) ÷ 2

    obj = formObjective(H, Z, Ω)

    function G(a, α, k)
        α_ = copy(α)
        α_[k] = a[1]
        return obj(α_)
    end

    α = copy(α0);
    
    for i = 1:n_iters, k = 1:(2*kmax)
        res = Optim.optimize(a -> G(a, α, k), amin, amax) # very slow and simple -- no gradient information

        α[k] = Optim.minimizer(res)[1]
    end
    
    return(α)
end


function formObjective(H, Z, Ω)
    ℓ       = maximum(Z)
    C       = evalCuts(H,Z,Ω)
    V, μ, S = evalSums(Z,H,ℓ,true);
    Ŝ       = aggregateSums(S,Ω)
    
    function objective(α)
        obj = 0.0
        for p in keys(Ŝ)
            Op   = Ω.ω(p, α)
            obj += get(C, p, 0)*log(Op) - Ŝ[p]*Op
        end
        return -convert(Float64, obj) # sign is for minimization
    end
    return objective
end