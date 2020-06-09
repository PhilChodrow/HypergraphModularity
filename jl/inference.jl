using StatsBase
using Combinatorics

include("HSBM.jl")
include("vol.jl")

function estimateΩ(H, Z; method="Piecewise Constant")
    @assert method == "Piecewise Constant"
    S = evalSums(Z,H)
    T = Dict()
    for k in keys(H.E), e in keys(H.E[k])
        z = Z[e]
        p = -sort(-collect(values(countmap(vec(z)))))
        n = values(countmap(vec(e)))
        m = values(countmap(p))
        T[p] = get(T,p,0) + H.E[k][e]
    end
    ω̂ = Dict()
    for p in keys(T)
       ω̂[p] = T[p]/S[p]
    end
    return ω̂
end