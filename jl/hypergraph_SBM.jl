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
    return(fill(S, rand(X)))
end

function sampleEdges(Z, ϑ, Ω, k_max, k_min=1)
    """
    run sampleEdge() for each possible distinct sequence of no more than k_max node labels, allowing repeats, and concatenate the results as a single array (edge list)
    The arrays Z and ϑ are required to be of the same length n. 
    """

    E = []
    n = length(Z) 
    for k = k_min:k_max
        T = with_replacement_combinations(1:n,k)
        for S = T
            e = sampleEdge(S, Z, ϑ, Ω)
            append!(E, e)
        end
    end
    return(E)
end

Z = [1 2 2 1 2]
ϑ = [4 3 5 2 1]

E = sampleEdges(Z, ϑ, z->plantedPartition(z,.001,.01), 3)
print(E)