function partitionize(a::Vector{<:Integer})
    """
    For a given integer vector a, return the partition corresponding to that
    vector. Useful for both counting corrections when sampling and computing
    likelihoods, and when computing partition-based values of Ω.

    This is the fastest version I could come up with.
    """
    a = sort(a)
    k = length(a)
    v = zero(a)

    v[1] = 1
    current = 1

    for i = 2:k
        if a[i] == a[i-1]
            v[current] += 1
        else
            current += 1
            v[current] = 1
        end
    end
    # v = v[v.>0]
    # return sort(v, rev = true)
    sort!(v, rev = true)
    return sortedremovezeros(v)
end

function sortedremovezeros(p::Vector{<:Integer})
    for i = 2:length(p)
        if p[i] == 0
            return p[1:i-1]
        end
    end
    return p
end

function counting_coefficient(z::Array{T, 1}) where {T<:Integer}
    p = partitionize(z)
    return Combinatorics.multinomial(p...)
end

function poisson_pdf(x::Integer, λ::Float64)
    exp(-λ)*λ^x/factorial(big(x))
end

function mutualInformation(Z, Ẑ, normalized = false)
    """
    Mutual information between two clusterings, optionally normalized
    Probably can be computed MUCH faster, but unlikely to be a bottleneck 
    in context. 
    """
    n = length(Z)
    
    p_XY = Dict()
    p_X = Dict()
    p_Y = Dict()

    for i = 1:length(Z)
        p_XY[(Z[i], Ẑ[i])] = get!(p_XY, (Z[i], Ẑ[i]), 0) + 1/n
        p_X[Z[i]]          = get!(p_X, Z[i], 0)          + 1/n
        p_Y[Ẑ[i]]          = get!(p_Y, Ẑ[i], 0)          + 1/n
    end
    
    I = 0
    for x in keys(p_X), y in keys(p_Y)
        try
            I += p_XY[x,y]*log(p_XY[x,y]/(p_X[x]*p_Y[y]))
        catch e
            nothing
        end
    end
    
    if normalized
        H_X, H_Y = 0, 0
        for x in keys(p_X)
            H_X -= log(p_X[x])*p_X[x]
        end
        for y in keys(p_Y)
            H_Y -= log(p_Y[y])*p_Y[y]
        end
        return (2*I)/(H_X + H_Y)
    end
    
    return I
end