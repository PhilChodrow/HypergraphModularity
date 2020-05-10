using Distributions
using Base 
using Combinatorics

include("omega.jl")

"""
The purpose of this module is to define a flexible stochastic blockmodel for hypergraphs. At this stage, all we are aiming to do is *sample* from the model give a partition, a group intensity function, and a degree vector. 
"""

function sampleEdge(S, Z, ϑ, Ω)
    """
    Sample a Poisson number of edges on a sequence of nodes S according to the hypergraph SBM law. 
    S: an array of node indices 
    Z: an integer array of length n. Each integer is a group label. 
    Ω: a group interaction function, such as plantedPartition
    Θ: a nonnegative float array of length n
    """
    z = Z[S]
    θ = ϑ[S]
    X = Poisson(prod(θ)*Ω(z))
    return(rand(X))
    # return(fill(S, rand(X)))
end

function sampleEdges(Z, ϑ, Ω; kmax=3, kmin=2)
    """
    run sampleEdge() for each possible distinct sequence of no more than k_max node labels, allowing repeats, and concatenate the results as a single array (edge list)
    The arrays Z and ϑ are required to be of the same length n. 
    Returns a Dict, keyed by edgesize. 
    Each value in this dict is itself a dict of edges, with counts, of the specified size. 
    """
    E = Dict{Integer, Dict}()
    n = length(Z)
    for k in kmin:kmax
        T = with_replacement_combinations(1:n, k)
        Ek = Dict{Array{Integer}, Integer}()
        for S in T
            X = sampleEdge(S, Z, ϑ, Ω)
            if X > 0
                Ek[S] = X
            end
        end
        E[k] = Ek
    end
    return(E)
end

function logLikelihood(E, Z, ϑ, Ω)
    """
    Given an edge list of the type returned by sampleEdges(), return the HSBM likelihood using labels Z, degree parameters ϑ, and group intensities Ω
    """
    n = length(Z)
    L = 0
    for k in keys(E)
        T = with_replacement_combinations(1:n, k) 
        Ek = E[k]   
        for S in T
            z = Z[S]
            θ = ϑ[S]
            X = Poisson(prod(θ)*Ω(z))
            
            m = get(Ek, S, 0)
            L += log(pdf(X, m))
        end
    end
    return(L)
end

# quick experiment

n = 100
Z = rand(1:5, n)
ϑ = ones(1,n) + rand(1,n)
Ω = z->plantedPartition(z,.0001,.001)

kmin = 1
kmax = 4

E = sampleEdges(Z, ϑ, Ω; kmax=kmax, kmin=kmin)

println("Generated data has ")
for k in kmin:kmax
    println("$(length(E[k])) edges of size $k")
end

ll = logLikelihood(E, Z, ϑ, Ω)
println("The likelihood of the true partition is $(round(ll, digits = 2))")

Z_wrong = rand(1:5, n)

ll = logLikelihood(E, Z_wrong, ϑ, Ω)

println("The likelihood of a random partition is $(round(ll, digits =2))")





