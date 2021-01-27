using CSV
using DataFrames

function parse_committees(chamber)
    DF = DataFrame!(CSV.File("$(chamber)_committees.csv"))
    
    sub = combine(groupby(DF, :new_id)) do sdf
               sdf[argmax(sdf.new_id), :]
               end

    labels = sub[:, "party"]

    hedges = Vector{Vector{Int64}}()
    for sub ∈ groupby(DF, [:session, :committee])
        hedge = sub[!, :new_id]
        sort!(hedge)
        push!(hedges, hedge)
    end

    names = ["$f $l" for (f, l) in zip(sub[!, "first"], sub[!, "last"])]

    dir = string("$(chamber)-committees")
    if !isdir(dir)
        mkdir(dir)
    end

    # hyperedge data
    open("$dir/hyperedges-$dir.txt", "w") do f
        for hedge in hedges
            write(f, join(hedge, ','))
            write(f, "\n")
        end
    end

    # names of the integer labels
    open("$dir/label-names-$dir.txt", "w") do f
        write(f, "Democrat\n")
        write(f, "Republican\n")        
    end
    
    # labels for each node
    open("$dir/node-labels-$dir.txt", "w") do f
        for party in labels
            write(f, "$(party)\n")
        end
    end

    # names of the nodes
    open("$dir/node-names-$dir.txt", "w") do f
        for name in names
            write(f, "$(name)\n")
        end
    end
end
