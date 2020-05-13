using Statistics

include("hypergraph_SBM.jl")
include("omega.jl")

# quick experiment

n = 20
Z = rand(1:5, n)
ϑ = ones(1,n) + rand(1,n)
μ = mean(ϑ)

# Ω = z->plantedPartition(z,1,10, k->(2*μ*k)^(-k))
Ω = z->groupSizePartition(z, k->(2*μ*k)^(-k))


kmin = 1
kmax = 4

E = sampleEdges(Z, ϑ, Ω; kmax=kmax, kmin=kmin)

println("Generated data has ")
for k in kmin:kmax
    println("$(length(E[k])) edges of size $k")
end

println("The degree sequence is ")
println(D(E))

ll = logLikelihood(E, Z, ϑ, Ω)
println("The log-likelihood of the true partition is $(round(ll, digits = 2))")

Z_wrong = rand(1:5, n)

ll = logLikelihood(E, Z_wrong, ϑ, Ω)

println("The log-likelihood of a random partition is $(round(ll, digits = 2))")







