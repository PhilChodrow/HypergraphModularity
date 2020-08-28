function estimateΩEmpirically(H, Z; min_val=0, aggregator=p->p)
    """
    aggregator groups partitions to new keys, upon which ω̂ is estimated (default: no aggregation)
    PC: not currently populating with a minimum value
    """

    T = Dict{Vector{Int64}, Float64}()

    # think this block is essentially tracking cuts
    # maybe Nate knows a faster way to do it?
    for k in keys(H.E), e in keys(H.E[k])
        z = Z[e]
        p = partitionize(z)
        T[p] = get(T,p,0) + H.E[k][e]
    end


    # Initialize

    r = maximum(k for k in keys(H.E)) # size of largest hyperedge
    ℓ = r
    S = evalSums(Z,H,ℓ, true)[3]
    T_agg = Dict{Any,Float64}()
    S_agg = Dict{Any,Float64}()
    for p in keys(S)
        agg_key = aggregator(p)
        T_agg[agg_key] = 0.0
        S_agg[agg_key] = 0.0
    end

    # Aggregate both sums over the aggregator keys
    for p in keys(S)
        agg_key = aggregator(p)
        T_agg[agg_key] += get(T, p, min_val)
        S_agg[agg_key] += S[p]
    end

    # Estimate (same for each partition that maps to the same key)
    ω̂ = Dict{Vector{Int64}, Float64}()
    for p in keys(S)
        agg_key = aggregator(p)
        if S_agg[agg_key] == 0
            ω̂[p] = min_val
        else
            ω̂[p] = T_agg[agg_key] / S_agg[agg_key]
        end
    end

    function Ω̂(x; α, mode="group")
        if mode == "group"
            return ω̂[partitionize(x)]
        elseif mode == "partition"
            return ω̂[x]
        end
    end
    return Ω̂
end

function coordinateAscent(H, Z, Ω, α0; n_iters = 10, amin = 0, amax = 10)
    
    kmax = length(α0) ÷ 2

    modularityObjective = formObjective(H, Z, Ω)

    function G(a, α, k)
        α_ = copy(α)
        α_[k] = a[1]
        return modularityObjective(α_)
    end

    α = copy(α0);
    
    for i = 1:n_iters, k = 1:(2*kmax)
        res = Optim.optimize(a -> G(a, α, k),  amin, amax) # very slow and simple -- no gradient information
        α[k] = Optim.minimizer(res)[1]
    end
    
    return(α)
end


function learnParameters(H, Z, Ω, α0; verbose = false, ftol_abs = 1e-6)
    """
    a more reliable optimization to learn parameters given a partition, using the COBYLA algorithm from NLopt. 
    returns both the optimal α and the objective value at that point. 
    """
    
    k = length(α0)
    
    obj_ = formObjective(H, Z, Ω)
    obj(x, grad) = obj_(x)
    
    opt = NLopt.Opt(:LN_COBYLA, k)
    opt.min_objective = obj
    opt.ftol_abs = ftol_abs
    
    (minf,minx,ret) = NLopt.optimize(opt, α0)
    numevals = opt.numevals;
    if verbose
        println("got $minf at $minx after $numevals iterations (returned $ret)")
    end
    return(minx, minf)
end