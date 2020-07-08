using StatsBase
using Combinatorics
include("vol.jl");

n = 50000

Z = rand(1:50, n)
D = rand(2:100, n)

r = 10 # maximum hyperedge size

@time V, μ, M = evalSums(Z, D, r;constants=true, bigInt=true);
