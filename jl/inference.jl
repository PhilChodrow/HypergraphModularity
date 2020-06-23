using StatsBase
using Combinatorics

include("HSBM.jl")
include("vol.jl")

function estimateΩ(H, Z; method="Piecewise Constant")
    @assert method == "Piecewise Constant"
    T = Dict{Vector{Int64}, Float64}()

    for k in keys(H.E), e in keys(H.E[k])
        z = Z[e]
        p = -sort(-collect(values(countmap(vec(z)))))
        n = values(countmap(vec(e)))
        m = values(countmap(p))
        T[p] = get(T,p,0) + H.E[k][e]
    end

    ω̂ = Dict{Vector{Int64},Float64}()
    S = evalSums(Z,H)[3]
    for p in keys(S)
       ω̂[p] = get(T, p, 0)/S[p] # doing a zero default value may not be the best choice statistically speaking. 
    end
    return ω̂
end
