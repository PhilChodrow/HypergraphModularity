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

function organize_by_size(Ω_dict::Dict{Any, Any})
    """
    Given a Dict() in which keys are partitions, returns a Vector{Dict{Array{Integer}, Any}}() in which the kth entry is the Dict() of partitions and values of fixed size. Can be used to obtain better performance. 
    """
    Om = Vector{Dict}()
    for k = 1:kmax
        omk = Dict(p => Ω_dict[p] for p in partitions(k))
        push!(Om, omk)
    end
    return Om
end

function buildΩ(Ω_dict::Dict{Any, Any}; by_size=true)
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

function ΩFromDictBySize(Ω_dict::Dict{Array{T,1}, Float64}) where {T<:Integer}
    """
    Create an partition-based intensity function Ω by passing a Dict() of partition-value pairs. 
    Behind the scenes, the intensity function returns values from a Vector of Dict{Array{T,1}, Float64}, where T<:Integer. 
    The kth element of this Vector is the set of partitions of size k. 
    Ω_dict::Dict{Array{T,1}, Float64} the intensity function dict
    return:: Ω, an intensity function which, when passed a key from Ω_dict, returns the corresponding value. 
    """
    Om = organize_by_size(Ω_dict)
    function Ω(p::Array{T,1}; mode = "group", k=0)::Float64 where {T<:Integer}
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

