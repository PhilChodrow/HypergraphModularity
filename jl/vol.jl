using StatsBase
using Combinatorics

include("HSBM.jl")


# Compute a complete set of sums using a fast recursive algorithm based on a simple combinatorial relationship between the required sums. 

# ------------------------------------------------------------------------------
# COMPUTATION OF SUMS FROM SCRATCH
# ------------------------------------------------------------------------------

function correctOvercounting(M::Dict, p::Array)
    """
    Utility function: second term in the recurrence in the notes. 
    M::Dict{Array{Integer, 1}, Integer}, the current Dict of constants being updated. 
    p::Array{Integer, 1}, the partition for which we are currently computing the new constant. 
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

function computeMoments(Z::Array, D::Array, r, ℓ=0)
    """
    Compute the volumes of groups in a hypergraph and the 1st through rth moments of those volumes. 
    Z::Array{Integer, 1}, the group labels. 
    D::Array{Integer, 1}, the degrees of nodes
    r::Integer, the size of the largest moment to compute (corresponds to kmax, the size of the largest hyperedge)
    ℓ::Integer, the number of groups. Defaults to the maximum integer value of Z
    return: V::Array{Integer, 1} the vector of volumes (one for each group) and μ::Array{Integer, 1} the vector of moments (1 through r)
    """
    if ℓ==0
        ℓ=maximum(Z)
    end

    V = [sum(D.*(Z .== j)) for j in 1:ℓ]
    μ = [sum(V.^i) for i = 1:r]
    return(V, μ)
end

function evalConstants(r)
    """
    Evaluate the combinatorial constants which appear as multipliers to the constants computed in evalSums. 
    r::Integer, the largest partition size to evaluate. Corresponds to kmax, the size of the largest hyperedge. 
    return C::Dict{Array{Integer, 1}, Integer}, a Dict() of constants, one for each partition p. 
    """
    C = Dict{Array{Integer, 1}, Integer}()
    for i = 1:r, j = 1:i, p in partitions(i, j)
        orderCorrection = prod([factorial(c) for c in values(countmap(p))])
        # orderCorrection = 1
        C[p] = multinomial(p...) ÷ orderCorrection
    end
    return(C)
end;

function evalSums(Z::Array, D::Array, r; constants=true, ℓ=0, bigInt=true)
    """
    Compute the sums required by the volume term of the modularity using a fast recursive method. 
    Z::Array{Integer, 1}, an Array of integer group labels
    D::Array{Integer, 1}, an Array of degrees
    r::Integer, the size of the largest partition to compute a constant for. Should correspond to kmax, the size of the largest hyperedge in an observed hypergraph. 
    constants::Bool, whether return the results multiplied by the combinatorial constants given by evalConstants(). Should be true in user-facing contexts and false when efficiently evaluating increments within algorithms. 
    ℓ::Integer the number of groups. Defaults to the maximum integer value of Z. 
    bigInt::Bool, whether to convert the degree sequence D into Array{bigInt,1} before processing. This is recommended for instances of even modest size, since the sums involved can become very large. Setting this parameter to false may produce faster runtimes and smaller memory allocations at the cost of silent errors due to integer overflow. 
    """

    if bigInt
        D = convert(Array{BigInt,1}, D)
        Z = convert(Array{BigInt,1}, Z)
    end

    V, μ = computeMoments(Z, D, r, ℓ)
    
    M = Dict{Array{Integer, 1}, Integer}()

    for i = 1:r, j = 1:i, p in partitions(i, j)
        M[p] = μ[p[end]]*get(M, p[1:(end-1)], 1) - correctOvercounting(M,p)
    end

    if constants
        C = evalConstants(r)
        for p in keys(M)
            M[p] *= C[p] 
        end
    end
    return V, μ, M
end

function evalSums(Z::Array, H::hypergraph, ℓ=0; bigInt=true)
    """
    An additional method for evalSums to operate on an object of type hypergraph. 
    H::hypergraph, the hypergraph on which to compute. 
    For other parameters and return values, see evalSums() above. 
    """
    r = maximum(keys(H.E))
    return evalSums(Z, H.D, r; ℓ=ℓ, bigInt=bigInt)
end

# ------------------------------------------------------------------------------
# INCREMENTS
# ------------------------------------------------------------------------------

function momentIncrements(V, μ, i, t, D, Z)
    """
    Compute the increment vectors in the volumes V and volume-moments μ associated with moving node i from its current group to group t. 
    V: an array of volumes associated with each group 
    μ: an array of moments of V
    i: the node to move
    t: the group into which node i will be moved
    D: the degree sequence of the hypergraph on which we compute
    Z: the current grouping vector
    return: ΔV and Δμ, the increments in V and μ
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

function increments(V, μ, M, i, t, D, Z)
    """
    compute increments in the sums M using the recursion formula from the main document. 
    V: an array of volumes associated with each group 
    μ: an array of moments of V
    M: the current Dict() of sums to be updated
    i: the node to move
    t: the group into which node i will be moved
    D: the degree sequence of the hypergraph on which we compute
    Z: the current grouping vector
    return ΔV, Δμ, ΔM: the increments in V, μ, and M
    """
    ΔV, Δμ = momentIncrements(V, μ, i, t, D, Z)

    # compute increments in M using recursion formula from notes
    ΔM = Dict{Array{Integer, 1}, Integer}()
    r = maximum([sum(p) for p in keys(M)])
    for i = 1:r, j = 1:i, p in partitions(i, j)
        ΔM[p] = Δμ[p[end]]*get(M, p[1:(end-1)], 1) + μ[p[end]]*get(ΔM, p[1:(end-1)], 0) + Δμ[p[end]]*get(ΔM, p[1:(end-1)], 0) - correctOvercounting(ΔM,p)
    end

    return(ΔV, Δμ, ΔM)
end

function addIncrements(V, μ, M, ΔV, Δμ, ΔM)
    """
    Add computed increments to the vector of volumes V, volume-moments μ, and sums M, obtaining updated versions. 
    V: an array of volumes associated with each group 
    μ: an array of moments of V
    M: the current Dict() of sums to be updated
    ΔV, Δμ, ΔM: increments in these vectors as computed by increments()
    return: V, μ, M, updated versions of each of the first three inputs. 
    """
    M̃ = Dict(p => M[p] + ΔM[p] for p in keys(M))
    return(V + ΔV, μ + Δμ, M̃)
end

# ------------------------------------------------------------------------------
# COMPLETE COMPUTATION OF SECOND (VOLUME) TERM IN MODULARITY
# ------------------------------------------------------------------------------

function second_term_eval(H::hypergraph, Z::Array{Int64, 1}, Ω; ℓ = 0, bigInt=true)
    """
    Compute the volume (second) term in hypergraph modularity. 
    H::hypergraph, the hypergraph on which to compute. 
    Z::Array{Integer, 1}, the vector of group labels. 
    Ω, the group interaction function. Needs to have a mode argument which, when set to value "partition", will cause evaluation on partition vectors rather than label vectors. 
    ℓ::Integer, the number of groups. Defaults to the maximum of the entries of Z. 
    bigInt::Bool, whether to convert the degrees to bigInts when performing summation. Will silently produce integer overflow errors on even modestly-sized instances. Only set to false if you KNOW that this won't be an issue. 
    """

    obj = 0
    kmax = maximum(keys(H.E))

    if ℓ == 0
        ℓ = maximum(Z)
    end

    V, μ, M = evalSums(Z, H, ℓ; bigInt=true)
    for p in keys(M)
        obj += Ω(p, mode = "partition")*M[p]
    end
    return obj
end