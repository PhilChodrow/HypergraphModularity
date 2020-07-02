using StatsBase
using Combinatorics

include("HSBM.jl")
include("vol.jl")

function estimateΩEmpirically(H, Z; min_val=0, aggregator=p->p)
    """
    Could organize by size later if we thought that would speed things up at all.

    aggregator groups partitions to new keys, upon which ω̂ is estimated (default: no aggregation)
    """
    
    T = Dict{Vector{Int64}, Float64}()

    for k in keys(H.E), e in keys(H.E[k])
        z = Z[e]
        p = partitionize(z)
        m = values(countmap(p))
        T[p] = get(T,p,0) + H.E[k][e]
    end

    # Initialize
    S = evalSums(Z,H)[3]
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
        ω̂[p] = T_agg[agg_key] / S_agg[agg_key]
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
