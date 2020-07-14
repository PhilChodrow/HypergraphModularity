## Generate a graph
using StatsBase
using Combinatorics

include("omega.jl")
include("HSBM.jl")
include("hypergraph_louvain.jl")
include("inference.jl");

n = 200
Z = rand(1:2, n)
ϑ = dropdims(ones(1,n) + rand(1,n), dims = 1)

# defining group intensity function Ω
μ = mean(ϑ)

ω(p,α) = (10 .*μ*sum(p))^(-sum(p))*prod(p.^α)^(1/(sum(p)*α))
α0 = 1

kmax = 3

Ω = buildΩ(ω, α0, kmax)
## Sample
H = sampleSBM(Z, ϑ, Ω; α=α0, kmax=kmax, kmin = 1)

# @time Z = SuperNodeLouvain(H,kmax,Ω;α=α0)

Juno.@profiler Ẑ = HyperLouvain(H,kmax,Ω;α=α0)

Juno.@profiler Ẑ = SuperNodeLouvain(H,kmax,Ω;α=α0)
