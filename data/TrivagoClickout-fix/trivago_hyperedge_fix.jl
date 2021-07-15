
f = readlines("hyperedges-TrivagoClickout-old.txt")


open("hyperedges-TrivagoClickout-old.txt", "w") do io
end

## parse
for line in f
    hyperedge = parse.(Int64,split(line,","))[1:end-1]
    open("hyperedges-TrivagoClickout-fix.txt", "a") do io
        for e = 1:length(hyperedge)-1
            write(io,"$(hyperedge[e]),")
        end
        write(io,"$(hyperedge[end])\n")
    end
end