using StatsBase

function partitionize2(z::Array{T, 1}) where {T<:Integer}
    """
    For a given integer vector z, return the partition corresponding to that vector. Useful for both counting corrections when sampling and computing likelihoods, and when computing partition-based values of Ω. 
    DEPRECATED: much slower than partitionize() below. 
    """
    p = collect(values(countmap(z)))
    sort!(p, rev=true)
end

function partitionize(a::Array{Int64,1})
    """
    For a given integer vector z, return the partition corresponding to that vector. Useful for both counting corrections when sampling and computing likelihoods, and when computing partition-based values of Ω. 
    """
    k = length(a)
    u = unique(a)
    d = Dict()
    lu = length(u)
    for j = 1:lu
        d[u[j]] = j
    end
    cnts = zeros(Int, lu)
    for i = 1:k
        @inbounds x = a[i]
        @inbounds cnts[d[x]] += 1
    end
    return sort(cnts,rev = true)
end


function counting_coefficient(z::Array{T, 1}) where {T<:Integer}
    p = partitionize(z)
    return multinomial(p...)
end

function poisson_pdf(x::Integer, λ::Float64)
    exp(-λ)*λ^x/factorial(x) 
end

    