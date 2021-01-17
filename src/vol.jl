# Compute a complete set of sums using a fast recursive algorithm based on a simple combinatorial relationship between the required sums.
# This code has been partially optimized for quick updates, which means that it is more complicated that it would need to be if we were just updating from scratch every time.

# Throughout these functions:
# - N::Dict{Vector{Int64},<:Integer} is a Dict mapping sorted partition vectors to their volume sums.
# - M::Dict{Vector{Int64},<:Integer} is a Dict mapping sorted partition vectors to their *uncorrected* sums. By *uncorrected*, we mean that these sums must be multipled by a set of combinatorial constants prior to including in the modularity function. Each entry of M differs from the corresponding entry of N by a combinatorial constant.
# - p::Vector{Int64} is a partition vector
# - D::Array{Int64, 1} is the degree sequence of the hypergraph, which does not change.
# - Z::Array{Int64, 1} is the group label vector, which does change.
# - ℓ::Int64 is the number of groups.

# ------------------------------------------------------------------------------
# COMPUTATION OF SUMS FROM SCRATCH
# ------------------------------------------------------------------------------

function correctOvercounting(M::Dict{Vector{Int64},<:Integer}, p::Vector{Int64})
    """
    Internal function: no reason to call this outside of evalSums().
    Computes the second term in the recurrence relation given in the paper.
    M::Dict{Vector{Int64},<:Integer}, a Dict mapping sorted partition vectors to their *uncorrected* volume-sums,
    p::Vector{Int64}, the partition for which we are calculating.
    returns S::bigInt, the second term in the recurrence relation.
    """
    pk = p[end]
    S = 0
    for i = 1:length(p)-1
        p_ = copy(p)[1:(end-1)]
        p_[i] += pk
        S += M[-sort(-p_)]
    end
    return(S)
end

function computeMoments(Z::Vector{<:Integer}, D::Vector{<:Integer}, r::Integer, ℓ::Integer=0)
    """
    Compute the vectors V of group volumes and μ of volume-moments.
    No reason to call this outside of evalSums()
    Z::Vector{Int64}, the vector of cluster labels.
    D::Vector{Int64}, the degree sequence
    r::Int64, the size of the largest hyperedge
    ℓ::Int64, the number of groups (defaults to maximum(Z))
    returns V::Array{Int64, 1} the vector of group volumes and μ::Vector{Int64} of volume-moments.
    NOTE: might want to implement bigInts in μ
    """
    if ℓ==0
        ℓ=maximum(Z)
    end

    V = [sum(D[Z .== j]) for j in 1:ℓ]
    μ = [sum(V.^i) for i = 1:r]
    return(V, μ)
end

function evalConstants(r::Int64; ℓ = r, bigInt::Bool=true)
    """
    Evaluate the combinatorial constants required to convert the "convenient" sums M to the required sums N required for modularity calculations.
    r::The size of the largest hyperedge.
    bigInt:Bool, whether to use bigInts for this computation, recommended unless the instance is VERY small.
    return: C::Dict{Array{Int64, 1}, bigInt} a Dict mapping a partition vector to its associated combinatorial constant.
    """

    C =
        if bigInt C = Dict{Vector{Int64},BigInt}()
        else      C = Dict{Vector{Int64},Int64}()
        end

    for i = 1:r, j = 1:minimum([i,ℓ]), p in Combinatorics.partitions(i, j)
        orderCorrection =
            if bigInt prod([factorial(big(c)) for c in values(StatsBase.countmap(p))])
            else      prod([factorial(c) for c in values(StatsBase.countmap(p))])
            end
        # orderCorrection = 1
        C[p] = Combinatorics.multinomial(big.(p)...) ÷ orderCorrection
    end
    return(C)
end;

function evalSums(Z::Vector{<:Integer}, D::Vector{Int64}, r::Integer; constants::Bool=true, ℓ::Integer=maximum(Z), bigInt::Bool=true)
    """
    Z::Array{Int64, 1} the vector of integer group labels
    D::Array{Int64, 1} the vector of degrees
    r::Int64, the largest hyperedge size to compute
    constants::Bool, whether to postmultiply by the combinatorial constants computed by evalConstants prior to returning the values. Choose true if computing modularity directly using the results, choose false if performing incremental updates a la Louvain
    ℓ::Int64, the number of groups (defaults to maximum(Z))
    bigInt::Bool, whether to use bigInt conversions for D and Z. Recommended unless the instance is VERY small.
    return (V::Array{bigInt,1}, μ::Array{bigInt,1}, M::Dict{Vector{Int64},bigInt}, with V and μ as described in computeMoments(), and M::Dict{Vector{Int64},bigInt} is the array of (optionally uncorrected) volume sums.
    """
    
    if bigInt
        D = convert(Array{BigInt,1}, D)
    end

    V, μ = computeMoments(Z, D, r, ℓ)

    M = bigInt ? Dict{Vector{Int64}, BigInt}() : Dict{Vector{Int64}, Int64}()
    
    ix = ℓ == 1 ? length(Z) : ℓ
    
    for i = 1:r, j = 1:minimum([i,ix]), p in Combinatorics.partitions(i, j)
        M[p] = μ[p[end]]*get(M, p[1:(end-1)], 1) - correctOvercounting(M,p)
    end

    if constants
        C = evalConstants(r; bigInt=bigInt, ℓ = ix)
        for p in keys(M)
            M[p] *= C[p]
        end
    end
    return V, μ, M
end

# Extra methods for evalSums for working with degree dicts and the custom hypergraph class.

function evalSums(Z::Vector{<:Integer}, H::hypergraph, ℓ::Integer=maximum(Z), bigInt::Bool=true, constants::Bool = true)
    r = maximum(keys(H.E))
    return evalSums(Z, H.D, r; constants = constants, ℓ = ℓ, bigInt = bigInt)
end

# ------------------------------------------------------------------------------
# INCREMENTS
# ------------------------------------------------------------------------------

function momentIncrements(V::Vector{<:Integer}, μ::Vector{<:Integer}, i::Int64, t::Int64, D::Vector{<:Integer}, Z::Vector{<:Integer})
    """
    Compute the update for the vectors V and μ associated with moving node i from its current group to group t.
    V::Array{Integer, 1}, the vector of group volumes.
    μ::Array{Integer, 1}, the vector of volume-moments.
    i::Int64, the node to move
    t::Int64, the proposed new group for node i
    D::Array{Integer, 1}, the degree vector
    Z::Array{Integer, 1}, the vector of group labels
    returns ΔV, Δμ, the required updates in V and μ
    """

    ΔV = zero(V)
    # update V
    s      = Z[i];
    ΔV[t] += D[i];
    ΔV[s] -= D[i];

    r = length(μ)
    Δμ = [(V[s] + ΔV[s])^j - V[s]^j + (V[t] + ΔV[t])^j - V[t]^j for j in 1:r]

    return ΔV, Δμ
end

function momentIncrements(V::Vector{<:Integer}, μ::Vector{<:Integer}, I_::Vector{<:Integer}, t::Integer, D::Vector{<:Integer}, Z::Vector{<:Integer})
    """
    Compute the update for the vectors V and μ associated with moving all the nodes listed in I from their group to group t.
    THIS CODE ASSUMES THAT THE NODES IN I ARE ALL IN THE SAME GROUP, AND WILL GIVE INCORRECT RESULTS OTHERWISE.
    V::Array{Integer, 1}, the vector of group volumes.
    μ::Array{Integer, 1}, the vector of volume-moments.
    I::Vector{Int64}, the list of multiple nodes to move
    t::Int64, the proposed new group for node i
    D::Array{Integer, 1}, the degree vector
    Z::Array{Integer, 1}, the vector of group labels
    returns ΔV, Δμ, the required updates in V and μ
    """

    ΔV = zero(V)
    # update V
    s  = Z[I_[1]];

    for i ∈ I_
        ΔV[t] += D[i];
        ΔV[s] -= D[i];
    end

    r = length(μ)
    Δμ = [(V[s] + ΔV[s])^j - V[s]^j + (V[t] + ΔV[t])^j - V[t]^j for j in 1:r]

    return ΔV, Δμ
end

function increments(V::Vector{<:Integer}, μ::Vector{<:Integer}, M::Dict{Vector{Int64},<:Integer}, I_::T, t::Integer, D::Vector{<:Integer}, Z::Vector{<:Integer}; bigInt::Bool=true) where T <: Union{<:Integer, Vector{<:Integer}}
    """
    Compute the update for the vectors V and μ, as well as the uncorrected sums M, associated with moving node(s) I from its current group to group t.
    THIS CODE ASSUMES THAT THE NODES IN I ARE ALL IN THE SAME GROUP, AND WILL GIVE INCORRECT RESULTS OTHERWISE.
    V::Array{Integer, 1}, the vector of group volumes.
    μ::Array{Integer, 1}, the vector of volume-moments.
    M::Dict{Array{Int64, 1}, bigInt}, the Dict of uncorrected volume sums.
    i::Int64, the node to move
    t::Int64, the proposed new group for node i
    D::Array{Integer, 1}, the degree vector
    Z::Array{Integer, 1}, the vector of group labels
    returns (ΔV, Δμ, ΔM), the required updates in V, μ, and M
    """

    ΔV, Δμ = momentIncrements(V, μ, I_, t, D, Z)

    # compute increments in M using recursion formula from notes
    if bigInt ΔM = Dict{Vector{Int64}, BigInt}()
    else      ΔM = Dict{Vector{Int64}, Int64}()
    end

    r = maximum([sum(p) for p in keys(M)])

    for i = 1:r, j = 1:i, p in Combinatorics.partitions(i, j)
        ΔM[p] = Δμ[p[end]]*get(M, p[1:(end-1)], 1) + μ[p[end]]*get(ΔM, p[1:(end-1)], 0) + Δμ[p[end]]*get(ΔM, p[1:(end-1)], 0) - correctOvercounting(ΔM,p)
    end
    return(ΔV, Δμ, ΔM)
end

function addIncrements(V::Vector{<:Integer}, μ::Vector{<:Integer}, M::Dict{Vector{Int64},<:Integer}, ΔV::Vector{<:Integer}, Δμ::Vector{<:Integer}, ΔM::Dict{Vector{Int64},<:Integer})
    """
    Add the increments (ΔV, Δμ, ΔM) to (V, μ, M), entrywise in the case of M.
    V::Array{Integer, 1}, the vector of group volumes.
    μ::Array{Integer, 1}, the vector of volume-moments.
    M::Dict{Array{Int64, 1}, bigInt}, the Dict of uncorrected volume sums.
    ΔV, Δμ, ΔM: increments in each of the above quantities returned by increments()
    """
    M̃ = Dict(p => M[p] + ΔM[p] for p in keys(M))
    return(V + ΔV, μ + Δμ, M̃)
end

# ------------------------------------------------------------------------------
# COMPLETE COMPUTATION OF SECOND (VOLUME) TERM IN MODULARITY
# ------------------------------------------------------------------------------

function aggregateSums(M, Ω::IntensityFunction)
    M̂ = Dict()
    for p in keys(M)
        p̂ = Ω.aggregator(p)
        M̂[p̂] = get(M̂, p̂, 0) + M[p]
    end
    return M̂
end



function second_term_eval(H::hypergraph, Z::Vector{<:Integer}, Ω::IntensityFunction; α, ℓ::Int64 = 0, bigInt::Bool=true)
    """
    Naive implementation, computes sums from scratch.
    H::hypergraph
    Z::Array{Int64, 1}, the group label vector.
    ℓ::Int64, number of groups
    Ω: group interation function (e.g., planted partition). Needs to have a mode argument which, when set to value "partition", will cause evaluation on partition vectors rather than label vectors.
    bigInt::Bool, whether to use bigInt conversions. Strongly recommended unless the instance is VERY small.
    """

    obj = 0.0

    if ℓ == 0
        ℓ = maximum(Z)
    end

    V, μ, M = evalSums(Z, H, ℓ, bigInt)
    M̂ = aggregateSums(M, Ω)
    for p̂ in keys(M̂)
        if M̂[p̂] >= 1
            obj += Ω.ω(p̂, α)*M̂[p̂]
        end
    end
    return obj
end
