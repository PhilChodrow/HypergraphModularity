using StatsBase

include("utils.jl")

"""
Throughout the docstrings, n gives the number of nodes. 
"""

export plantedPartition

function harmonicMean(p::Array{T, 1}) where {T<:Integer}
    """
    return the harmonic mean of array p; used for some toy interaction functions. 
    """
    k = sum(p)
    prod(p)^(1/k)
end

function organize_by_size(Ω_dict::Dict{Vector{Int64}, Float64})
    """
    Given a Dict() in which keys are partitions, returns a Vector{Dict{Array{Integer}, Any}}() in which the kth entry is the Dict() of partitions and values of fixed size. Can be used to obtain better performance. 
    """
    Om = Dict{Int64, Dict{Vector{Int64}, Float64}}()
    for p in keys(Ω_dict)
        k = sum(p)
        if !haskey(Om,k)
            Om[k] = Dict()
        end
        Om[k][p] = Ω_dict[p]
    end
    return Om        
    # end
    # for k = 1:kmax
    #     omk = Dict(p => Ω_dict[p] for p in partitions(k))
    #     push!(Om, omk)
    # end
    # return Om
end

function buildΩ(Ω_dict::Dict{Array{T,1}, Float64}; by_size=true) where {T<:Integer}
    """
    D: a Dict() in which the keys are partition vectors. 
    returns: Ω, an interaction function which when evaluated on p returns D[p] if p is a partition vector, mode = "partition" or D[partitionize(p)] if p is a list of group labels (mode = "group"). 
    Only recommended for use with generalized partition-based models. 
    Could be generalized but hasn't been yet. 
    """
    # if by_size is true, the output function will optionally allow the specification of the edge size k, saving the need to recompute it from each partition size. 
    if by_size
        return ΩFromDictBySize(Ω_dict)
    else
        return ΩFromDict(Ω_dict)
    end
end


function ΩFromDict(Ω_dict::Dict{Array{T,1}, Float64}) where {T<:Integer}
    """
    Create an partition-based intensity function Ω by passing a Dict() of partition-value pairs. 
    Ω_dict::Dict{Array{T,1}, Float64} the intensity function dict
    return:: Ω, an intensity function which, when passed a key from Ω_dict, returns the corresponding value. 
    """
    function Ω(p::Array{T,1}; mode="group", k=0)::Float64 where {T<:Integer}
        if mode == "group"
            return Ω_dict[partitionize(p)]
        elseif mode == "partition"
            return Ω_dict[p]
        end
    end
    return Ω
end

function ΩFromDictBySize(Ω_dict::Dict{Vector{Int64}, Float64}) 
    """
    Create an partition-based intensity function Ω by passing a Dict() of partition-value pairs. 
    Behind the scenes, the intensity function returns values from a Vector of Dict{Array{T,1}, Float64}, where T<:Integer. 
    The kth element of this Vector is the set of partitions of size k. 
    Ω_dict::Dict{Array{T,1}, Float64} the intensity function dict
    return:: Ω, an intensity function which, when passed a key from Ω_dict, returns the corresponding value. 
    """
    Om = organize_by_size(Ω_dict)
    function Ω(p::Vector{<:Integer}; mode = "group", k=0)::Float64 
        if mode == "group"
            if k == 0
                k = length(p)
            end
            return Om[k][partitionize(p)]
        elseif mode == "partition"
            if k == 0
                k = sum(p)
            end
            return Om[k][p]
        end
    end
    return Ω
end

function buildΩ(f, α0, kmax)
    """
    f : function, takes in *partitions* and returns real numbers. May need to check on whether Float64 is big enough
    α0: parameters passed to f
    kmax: the size of the largest collection of nodes for which Ω should be able to compute
    """
    # memoization
    ᾱ = α0
    Om = Dict(k => Dict(p => f(p,ᾱ) for p in partitions(k)) for k = 1:kmax)
    function Ω(x; α, mode="group")
        
        # check if the parameter has been updated and re-memoize if so
        if !(α ≈ ᾱ)
#             println("rememoizing")
            ᾱ = α
            Om = Dict(k => Dict(p => f(p,ᾱ) for p in partitions(k)) for k = 1:kmax)
        end
        
        # otherwise, return required value from Om
        if mode == "group"
            return Om[length(x)][partitionize(x)]
        elseif mode == "partition"
            return Om[sum(x)][x]
        end
    end
    return Ω
end


