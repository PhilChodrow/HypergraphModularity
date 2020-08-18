"""
Throughout the docstrings, n gives the number of nodes.
"""

mutable struct IntensityFunction
    """
    range should be a vector of all valid inputs into Ω
    Range of P should be range
    """
    ω
    P
    range
    aggregator
end

function partitionsUpTo(kmax)
    return [p for k = 1:kmax for p in Combinatorics.partitions(k) ]
end

function partitionIntensityFunction(ω, kmax)
    range      = partitionsUpTo(kmax)
    P          = partitionize
    aggregator = identity
    return IntensityFunction(ω, P, range, aggregator)
end

function allOrNothingIntensityFunction(ω, kmax)
    range      = [(1.0*x, y) for x = 0:1 for y = 1:kmax]
    P          = z->(all(z[1] .== z), length(z))
    aggregator = p->(length(p) == 1, sum(p))
    return IntensityFunction(ω, P, range, aggregator)
end

function sumOfExteriorDegreesIntensityFunction(ω, kmax)
    range = [(1.0*x, y) for y = 1:kmax for x = 0:y]

    function P(z)
        len = length(z)
        return (length(unique!(z)), len)
    end

    aggregator = p->(length(p), sum(p))
    return IntensityFunction(ω, P, range, aggregator)
end







# function buildΩ(f, α0, kmax)
#     """
#     f : function, takes in *partitions* and returns real numbers. May need to check on whether Float64 is big enough
#     α0: parameters passed to f
#     kmax: the size of the largest collection of nodes for which Ω should be able to compute
#     """
#     # memoization
#     ᾱ = α0
#     Om = Dict(k => Dict(p => f(p,ᾱ) for p in Combinatorics.partitions(k)) for k = 1:kmax)
#     function Ω(x; α, mode="group")
#
#         # check if the parameter has been updated and re-memoize if so
#         if !(α ≈ ᾱ)
#             ᾱ = α
#             Om = Dict(k => Dict(p => f(p,ᾱ) for p in Combinatorics.partitions(k)) for k = 1:kmax)
#         end
#
#         # otherwise, return required value from Om
#         if mode == "group"
#             return Om[length(x)][partitionize(x)]
#         elseif mode == "partition"
#             return Om[sum(x)][x]
#         end
#     end
#     return Ω
# end
