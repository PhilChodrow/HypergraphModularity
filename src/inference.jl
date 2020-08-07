using StatsBase
using Combinatorics

include("HSBM.jl")
include("vol.jl")

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

function estimateParameters(H, Z, Ω, α0)
    """
    Currently assumes that the parameter vector has 2*kmax entries, and that the first and second halves of the parameter vector mean meaningfully different things.
    VERY JANKY ATM.
    """

    ℓ = maximum(k for k in keys(H.E)) # size of largest hyperedge

    C       = evalCuts(Z,H)
    V, μ, S = evalSums(Z,H,ℓ,true);

    α = copy(α0)
    kmax = length(α)÷2

    res = 0

    function objective(α)
        obj = 0.0
        for p in keys(S)
            Op   = Ω(p; α=α, mode="partition")
            obj += get(C, p, 0)*log(Op) - S[p]*Op
        end
        return -Float64(obj, RoundDown) # sign is for minimization
    end

    function objective(α, a, k)
        α_ = copy(α)
        α_[k] = a[1]
        return objective(α_)
    end

    for i = 1:50
#         println(-Optim.minimum(res))
        # optimization in γ
        for k = (kmax+1):(2*kmax)
            res = optimize(a -> objective(α, a, k), -100.0, 100.0) # very slow and simple -- no gradient information
            α[k] = Optim.minimizer(res)[1]

        end
        # optimizationin β
        for k = 1:kmax
            res = optimize(a -> objective(α, a, k), -100, 100)
            α[k] = Optim.minimizer(res)[1]
        end
    end

    return α, -Optim.minimum(res)
end
