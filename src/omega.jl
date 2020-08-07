using StatsBase

include("utils.jl")

"""
Throughout the docstrings, n gives the number of nodes. 
"""

export plantedPartition

function buildΩ(f, α0, kmax)
    """
    f : function, takes in *partitions* and returns real numbers. May need to check on whether Float64 is big enough
    α0: parameters passed to f
    kmax: the size of the largest collection of nodes for which Ω should be able to compute
    """
    # memoization
    ᾱ = α0
    Om = Dict(k => Dict(p => f(p,ᾱ) for p in partitions(k)) for k = 1:kmax)
    function Ω(x; α, mode="group")
        
        # check if the parameter has been updated and re-memoize if so
        if !(α ≈ ᾱ)
            ᾱ = α
            Om = Dict(k => Dict(p => f(p,ᾱ) for p in partitions(k)) for k = 1:kmax)
        end
        
        # otherwise, return required value from Om
        if mode == "group"
            return Om[length(x)][partitionize(x)]
        elseif mode == "partition"
            return Om[sum(x)][x]
        end
    end
    return Ω
end


