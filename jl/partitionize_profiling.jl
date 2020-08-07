include("utils.jl")

nreps = 100000

using DataStructures

Z = [3, 2, 2,2, 1, 3, 4, 3, 4, 4, 4, 5, 4, 5, 4, 5, 1]
# Z = [20, 4, 1, 1]
Z_ = convert(Array{Int16, 1}, Z)


function partitionize_phil(a::Vector{<:Integer})
    a = sort(a)
    k = length(a)
    v = zero(a)

    v[1] = 1
    current = 1

    for i = 2:k
        if a[i] == a[i-1]
            v[current] += 1
        else
            current += 1
            v[current] = 1
        end
    end
    v = v[v.>0]
    return sort(v, rev = true)
end


function partitionize_nate(a::Array{<:Integer,1})
    k = length(a)
    d = Dict{Int64,Int64}()
    d[a[1]] = 1     # dictionary from cluster index to id in vector v
    next = 2
    v = [1]         # number of times the cluster index shows up
    for i = 2:k
        ind = get(d, a[i], next)
        if ind == next
            push!(v,1)
            next += 1
        else
            v[ind] += 1
        end
    end
    return sort(v,rev = true)
end

@time for i = 1:nreps partitionize(Z) end
@time for i = 1:nreps partitionize(Z_) end

@time for i = 1:nreps partitionize_nate(Z) end
@time for i = 1:nreps partitionize_phil(Z) end

@time for i = 1:nreps partitionize_nate(Z_) end
@time for i = 1:nreps partitionize_phil(Z_) end






Juno.@profiler for i = 1:nreps partitionize(Z) end
Juno.@profiler for i = 1:nreps partitionize(Z_) end

Juno.@profiler for i = 1:nreps partitionize_nate(Z) end
Juno.@profiler for i = 1:nreps partitionize_phil(Z) end
