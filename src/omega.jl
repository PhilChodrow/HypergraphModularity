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
    range      = [(x, y) for x = 0:1 for y = 1:kmax]
    P          = z->(all(x->x==z[1],z), length(z))
    aggregator = p->(length(p) == 1, sum(p))
    return IntensityFunction(ω, P, range, aggregator)
end

function sumOfExteriorDegreesIntensityFunction(ω, kmax)
    range = [(1*x, y) for y = 1:kmax for x = 0:y]

    function P(z)
        len = length(z)
        return (length(unique!(z)), len)
    end

    aggregator = p->(length(p), sum(p))
    return IntensityFunction(ω, P, range, aggregator)
end