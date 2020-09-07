function read_hypergraph_label_names(dataname::String)
    names = String[]
    open("data/$dataname/label-names-$dataname.txt") do f
        for line in eachline(f)
            push!(names, line)
        end
    end
    return names
end

function read_hypergraph_labels(dataname::String)
    labels = Int64[]
    open("data/$dataname/node-labels-$dataname.txt") do f
        for line in eachline(f)
            push!(labels, parse(Int64, line))
        end
    end
    return labels
end

function read_hypergraph_edges(dataname::String, maxsize::Int64=25)
    E = Dict{Integer, Dict}()
    open("data/$dataname/hyperedges-$dataname.txt") do f
        for line in eachline(f)
            edge = [parse(Int64, v) for v in split(line, ',')]
            sort!(edge)
            if length(edge) > maxsize; continue; end
            sz = length(edge)
            if !haskey(E, sz)
                E[sz] = Dict{}()
            end
            E[sz][edge] = 1
        end
    end
    return E
end

function read_hypergraph_data(dataname::String, maxsize::Int64=25)
    labels = read_hypergraph_labels(dataname)
    E = read_hypergraph_edges(dataname, maxsize)
    
    n = length(labels)
    D = zeros(Int64, n)
    for (sz, edges) in E
        for (e, _) in edges
            D[e] .+= 1
        end
    end

    N = 1:n
    return hypergraph(N, E, D), labels
end
;
