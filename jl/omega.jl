using StatsBase

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

function groupSizePartition(z, fk=k->1)
    k = length(z)
    a = countmap(vec(z))
    a = sort(collect(values(a)))
    fk(k)*prod(a)^(1/k) 
end
