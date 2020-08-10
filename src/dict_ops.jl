import Base: *, +, -, /

function *(d1::Dict, d2::Dict)
    """
    entrywise multiplication. Will error if any keys are present in d1 and not in d2
    """
    d3 = Dict(p => 0 for p in keys(d1))
    for p in keys(d1)
        d3[p] = d1[p]*d2[p]
    end
    return d3
end

function +(d1::Dict, d2::Dict)
    """
    entrywise addition. Will error if any keys are present in d1 and not in d2
    """
    d3 = Dict(p => 0 for p in keys(d1))
    for p in keys(d1)
        d3[p] = d1[p]+d2[p]
    end
    return d3
end

function -(d1::Dict, d2::Dict)
    """
    entrywise subtraction. Will error if any keys are present in d1 and not in d2
    """
    d3 = Dict(p => 0 for p in keys(d1))
    for p in keys(d1)
        d3[p] = d1[p]-d2[p]
    end
    return d3
end


function /(d1::Dict, d2::Dict)
    """
    entrywise division. Will error if any keys are present in d1 and not in d2
    """
    d3 = Dict(p => 0 for p in keys(d1))
    for p in keys(d1)
        d3[p] = d1[p]/d2[p]
    end
    return d3
end
