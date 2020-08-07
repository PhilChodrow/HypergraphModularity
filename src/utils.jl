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
    v = v[v.>0]
    return sort(v, rev = true)
end

function counting_coefficient(z::Array{T, 1}) where {T<:Integer}
    p = partitionize(z)
    return Combinatorics.multinomial(p...)
end

function poisson_pdf(x::Integer, λ::Float64)
    exp(-λ)*λ^x/factorial(x)
end
