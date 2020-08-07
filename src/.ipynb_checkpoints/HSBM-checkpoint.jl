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
    E::Dict{Integer, Dict} 
    D::Dict{Integer, Integer} = Dict{Integer, Integer}()
end

function sampleEdge(S::Array{Int64,1}, Z::Array{Int64,1}, ϑ::Array{Float64,1}, Ω::Any)
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

function sampleEdges(Z::Array{Int64,1}, ϑ::Array{Float64,1}, Ω::Any; kmax::Integer=3, kmin::Integer=2)
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
        Ek = Dict{Array{Int64}, Integer}()
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

function D(E::Dict{Integer, Dict})
    """
    Compute the degree sequence of an edge list. 
    """
    d = Dict{Integer, Integer}()
    for k in keys(E)
        Ek = E[k]
        for e in keys(Ek)
            for i in e
                d[i] = get(d, i, 0) + 1
            end
        end
    end
    return(d)
end

function computeDegrees!(H::hypergraph)
    """
    Compute the degree sequence of a hypergraph and store it as a field of the hypergraph. 
    """
    H.D = D(H.E)
end

function sampleSBM(args...;kwargs...)
    """
    Sample a hypergraph with specified parameters and return it with its degree sequence pre-computed. This is the primary user-facing function for sampling tasks. 
    """
    E = sampleEdges(args...;kwargs...)
    H = hypergraph(E = E)

    computeDegrees!(H)
    return(H)
end

function logLikelihood(H::hypergraph, Z::Array{Int64,1}, ϑ::Array{Float64,1}, Ω::Any)
    """
    Given a hypergraph, return the HSBM likelihood using labels Z, degree parameters ϑ, and group intensities Ω.
    NOTE: this is a VERY slow function that should be spead up by orders of magnitude when Ω falls into important special cases cases. 
    """
    n = length(Z)
    L = 0
    for k in keys(H.E)  
        T = with_replacement_combinations(1:n, k) 
        Ek = H.E[k]   
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

function countEdges(H::hypergraph)
    """
    count the number of edges in H
    """
    sum([length(H.E[k]) for k in keys(H.E)])
end