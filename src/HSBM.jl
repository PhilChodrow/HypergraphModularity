"""
The purpose of this module is to define a flexible stochastic blockmodel for hypergraphs. At this stage, all we are aiming to do is *sample* from the model give a partition, a group intensity function, and a degree vector.
"""

Parameters.@with_kw mutable struct hypergraph
    """
    A very simple hypergraph composite type, designed to hold a node list N, an edge list E, a degree sequence D,
    """

    N::Vector{Int64}
    E::Dict{Int64, Dict}
    D::Array{Int64, 1} = Array{Int64, 1}()

end

function sampleEdge(S::Vector{Int64}, Z::Vector{Int64}, ϑ::Vector{Float64}, Ω::IntensityFunction; α::Any)
    """
    Sample a Poisson number of edges on a sequence of nodes S according to the hypergraph SBM law.

    # Arguments

    S::Vector{Int64},  an array of node indices
    Z::Vector{Int64}, an integer array of length n. Each integer is a group label.
    Ω::IntensityFunction, an IntensityFunction governing interaction rates between nodes of different groups. 
    ϑ::Vector{Float64},  a nonnegative float array of length n

    # Returns

    x::Int64, an integer number of edge counts sampled according to a Poisson random variable with specified formula. 

    """
    z = Z[S]
    θ = ϑ[S]

    # combinatorial factor associated with repeated indices. Equivalent to sampling a separate Poisson for each of the c distinct permutations of the node labels S
    c = counting_coefficient(S)
    X = Distributions.Poisson(prod(θ)*Ω.ω(Ω.P(z), α)*c)
    return(rand(X))
end

function sampleEdges(Z::Vector{Int64}, ϑ::Vector{Float64}, Ω::IntensityFunction; α::Vector{Float64}, kmax::Int64=3, kmin::Int64=2)
    """

    Sample edges of specified sizes according to the degree-corrected hypergraph configuration model (DCHSBM). Runs sampleEdge() for each possible distinct sequence of nodes of size between kmin and kmax, inclusive. 

    # Arguments

    Z::Vector{Int64}, a vector of group labels
    ϑ::Vector{Float64}, a vector of degree parameters
    Ω::IntensityFunction, an intensity function governing the rate of interaction between nodes of various groups.
    α::Vector{Float64}, a parameter passed to Ω
    kmax::Int64, the size of the largest hyperedge to form
    kmin::Int64, the size of the smallest hyperedge to form

    # Returns

    E::Dict{Int64,Dict{Vector{Int64}, Int64}}, a dictionary keyed by edge sizes. E[k] is itself a dictionary whose keys are ordered sets of nodes and whose values are the number of edges on those node sets. 
    N::Vector{Int64}, the list of nodes. 
    """
    E = Dict{Int64, Dict}() # initialize empty edge list
    n = length(Z)           # number of nodes inferred from group labels

    for k in kmin:kmax
        T = Combinatorics.with_replacement_combinations(1:n, k) # all distinct combos of k nodes
        Ek = Dict{Vector{Int64}, Int64}()         # list of edges of size k
        for S in T
            X = sampleEdge(S, Z, ϑ, Ω; α=α)
            if X > 0                              # only store combos with at least one edge
                Ek[S] = X
            end
        end
        E[k] = Ek
    end
    N = 1:n
    return(E, N)
end

function sampleEdges(Z::Dict, ϑ::Dict, Ω::IntensityFunction; α::Any, kmax::Int64=3, kmin::Int64=2)
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

    # Arguments

    args and kwargs: See documentation for sampleEdges()

    # Returns

    H::hypergraph, a hypergraph sampled according to the DCHSBM with specified args and kwargs. 
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


Base.copy(H::hypergraph) = hypergraph(H.N, H.E, H.D)