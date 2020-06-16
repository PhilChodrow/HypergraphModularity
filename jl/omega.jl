using StatsBase

include("utils.jl")

"""
Throughout the docstrings, n gives the number of nodes. 
"""

export plantedPartition, groupSizePartition

# test if all elements are equal, from 
# https://stackoverflow.com/questions/47564825/check-if-all-the-elements-of-a-julia-array-are-equal/47578613

δ(x) = all(y->y==x[1],x)

function plantedPartition(z, ω0, ω1, fk=k->1)
    """
    z: an array of nonnegative integers
    ω0: nonnegative float
    ω1: nonnegative float
    RETURN: ω1 if all entries of z are equal and ω0 otherwise
    """
    k = length(z)
    fk(k)*(ω0 + (ω1 - ω0)*δ(z))
end

function harmonicMean(p)
    k = sum(p)
    prod(p)^(1/k)
end

function ΩFromDict(D)
    """
    D: a Dict() in which the keys are partition vectors. 
    returns: Ω, an interaction function which when evaluated on p returns D[p] if p is a partition vector, mode = "partition" or D[partitionize(p)] if p is a list of group labels (mode = "group"). 
    Only recommended for use with generalized partition-based models. 
    Could be generalized but hasn't been yet. 
    """
    function Ω(p; mode="group")
        if mode == "group"
            return D[partitionize(p)]
        elseif mode == "partition"
            return D[p]
        end
    end
    return Ω
end


