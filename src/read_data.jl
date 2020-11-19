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

function read_hypergraph_data(dataname::String, maxsize::Int64=25, return_labels = true)
    
    E = read_hypergraph_edges(dataname, maxsize)
    
    n = maximum([maximum(e) for k in keys(E) for e in keys(E[k])])
    D = zeros(Int64, n)
    for (sz, edges) in E
        for (e, _) in edges
            D[e] .+= 1
        end
    end
    
    maxedges = maximum(keys(E))
    for k in 1:maxedges
        if !haskey(E, k)
            E[k] = Dict{}()
        end
    end
    
    N = 1:n
    
    if return_labels
        labels = read_hypergraph_labels(dataname)
        return hypergraph(N, E, D), labels
    end
    
    return hypergraph(N, E, D)
end
;

function hyperedges(dataname::String)
    nverts = Int64[]
    open("$(dataname)-nverts.txt") do f
        for line in eachline(f); push!(nverts, parse(Int64, line)); end
    end
    simps = Int64[]
    open("$(dataname)-simplices.txt") do f
        for line in eachline(f); push!(simps, parse(Int64, line)); end
    end
    edges = Set{Set{Int64}}()
    ind = 0
    for nvert in nverts
        edge = Set{Int64}(simps[(ind + 1):(ind + nvert)])
        push!(edges, edge)
        ind += nvert
    end
    open("hyperedges-$(dataname).txt", "w") do f
        for edge in edges
            vedge = collect(edge)
            sort!(vedge)
            s = join(vedge, ",")
            write(f, "$s\n")
        end
    end
end

