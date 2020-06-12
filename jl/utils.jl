using StatsBase

function partitionize(z)
    """
    For a given integer vector z, return the partition corresponding to that vector. Useful for both counting corrections when sampling and computing likelihoods, and when computing partition-based values of Ω. 
    It may be useful to replace this with cvec_2_pvec(), which currently lives in cut.jl
    """
    p = countmap(vec(z))
    p = -sort(-collect(values(p)))
end

function counting_coefficient(z)
    p = partitionize(z)
    return multinomial(p...)
end

function poisson_pdf(x, λ)
    exp(-λ)*λ^x/factorial(x) 
end

    