using Distributions
using Base 
using Combinatorics
using Parameters

include("omega.jl")


"""
The purpose of this module is to define a flexible stochastic blockmodel for hypergraphs. At this stage, all we are aiming to do is *sample* from the model give a partition, a group intensity function, and a degree vector. 
"""

@with_kw mutable struct hypergraph
    """
    A very simple hypergraph composite type, designed to hold an edge list E, a degree sequence. 
    """

    N::Vector{Int64}
    E::Dict{Int64, Dict} 
    D::Array{Int64, 1} = Array{Int64, 1}()

end

function sampleEdge(S::Vector{Int64}, Z::Vector{Int64}, ϑ::Vector{Float64}, Ω::Any; α::Any)
    """
    Sample a Poisson number of edges on a sequence of nodes S according to the hypergraph SBM law. 
    S: an array of node indices 
    Z: an integer array of length n. Each integer is a group label. 
    Ω: a group interaction function, such as plantedPartition
    Θ: a nonnegative float array of length n
    C: combinatorial factor
    """
    z = Z[S]
    θ = ϑ[S]
    
    # combinatorial factor associated with repeated indices. Equivalent to sampling a separate Poisson for each of the c distinct permutations of the node labels S
    c = counting_coefficient(S)    
    X = Poisson(prod(θ)*Ω(z; α=α, mode="group")*c)
    return(rand(X))
end

function sampleEdges(Z::Vector{Int64}, ϑ::Vector{Float64}, Ω::Any; α::Any, kmax::Int64=3, kmin::Int64=2)
    """
    run sampleEdge() for each possible distinct sequence of no more than k_max node labels, allowing repeats, and concatenate the results as a single array (edge list)
    The arrays Z and ϑ are required to be of the same length n. 
    Returns a Dict, keyed by edgesize. 
    Each value in this dict is itself a dict of edges, with counts, of the specified size. 
    """
    E = Dict{Int64, Dict}()
    n = length(Z)
    for k in kmin:kmax
        T = with_replacement_combinations(1:n, k)
        Ek = Dict{Vector{Int64}, Int64}()
        for S in T
            X = sampleEdge(S, Z, ϑ, Ω; α=α)
            if X > 0
                Ek[S] = X
            end
        end
        E[k] = Ek
    end
    N = 1:n
    return(E, N)
end

function sampleEdges(Z::Dict, ϑ::Dict, Ω::Any; α::Any, kmax::Int64=3, kmin::Int64=2)
    """
    An alternative method when Z and ϑ are specified as dicts keyed by node rather than as arrays. 
    """
    Z = [Z[i] for i in 1:length(Z)]
    ϑ = [ϑ[i] for i in 1:length(ϑ)]
    sampleEdges(Z, ϑ, Ω; α=α, kmax=kmax, kmin=kmin)
end

function computeDegrees(E::Dict{Int64, Dict}, N::Vector{Int64})
    """
    Compute the degree sequence of an edge list. 
    """

    d = zeros(length(N))

    for k in keys(E)
        Ek = E[k]
        for e in keys(Ek)
            for i in e
                d[i] += 1
            end
        end
    end
    return(d)

    # d = Dict(i => 0 for i in N)

    # for k in keys(E)
    #     Ek = E[k]
    #     for e in keys(Ek)   
    #         for i in e
    #             d[i] = get(d, i, 0) + 1
    #         end
    #     end
    # end
    # for i = 1:maximum(keys(d)) # fill in zeros 
    #     d[i] = get(d, i, 0)
    # end
    # return([d[i] for i = 1:length(d)])
end

function computeDegrees(H::hypergraph)
    return computeDegrees(H.E, H.N)
end

function computeDegrees!(H::hypergraph)
    """
    Compute the degree sequence of a hypergraph and store it as a field of the hypergraph. 
    """
    H.D = computeDegrees(H)
end

function sampleSBM(args...;kwargs...)
    """
    Sample a hypergraph with specified parameters and return it with its degree sequence pre-computed. This is the primary user-facing function for sampling tasks. 
    """
    E, N = sampleEdges(args...;kwargs...)
    H = hypergraph(E = E, N = N)

    computeDegrees!(H)
    return(H)
end

function countEdges(H::hypergraph)
    """
    count the number of edges in H
    """
    sum([length(H.E[k]) for k in keys(H.E)])
end
