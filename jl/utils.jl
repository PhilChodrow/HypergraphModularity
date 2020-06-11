using StatsBase

function partitionize(z)
    p = countmap(vec(z))
    p = -sort(-collect(values(p)))
end