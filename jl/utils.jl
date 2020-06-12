using StatsBase

function partitionize(z)
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

    