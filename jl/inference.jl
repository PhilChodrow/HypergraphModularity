using StatsBase
using Combinatorics

include("HSBM.jl")
include("vol.jl")

function estimateΩEmpirically(H, Z; min_val=0)
    """
    Could organize by size later if we thought that would speed things up at all.
    """
    
    T = Dict{Vector{Int64}, Float64}()

    for k in keys(H.E), e in keys(H.E[k])
        z = Z[e]
        p = partitionize(z)
        m = values(countmap(p))
        T[p] = get(T,p,0) + H.E[k][e]
    end

    ω̂ = Dict{Vector{Int64},Float64}()
    S = evalSums(Z,H)[3]
    for p in keys(S)
       ω̂[p] = get(T, p, min_val)/S[p] 
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
