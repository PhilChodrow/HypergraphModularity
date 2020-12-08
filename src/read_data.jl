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

function read_hypergraph_edges(dataname::String, maxsize::Int64=25, minsize::Int64=2)
    E = Dict{Integer, Dict}()
    open("data/$dataname/hyperedges-$dataname.txt") do f
        for line in eachline(f)
            edge = [parse(Int64, v) for v in split(line, ',')]
            sort!(edge)
            if minsize <= length(edge) <= maxsize
                sz = length(edge)
                if !haskey(E, sz)
                    E[sz] = Dict{}()
                end
                E[sz][edge] = 1
            end
        end
    end
    return E
end


function read_hypergraph_data(dataname::String, maxsize::Int64=25, minsize::Int64=2, return_labels = true)
    
    E = read_hypergraph_edges(dataname, maxsize, minsize)
    
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

function read_stats_data(kmax=0)
    path = "data/SCC2016-with-abs/SCC2016/Data/authorPaperBiadj.txt"
    A = DelimitedFiles.readdlm(path, '\t', Int, '\n');
    
    n = size(A)[1]
    N = 1:n
    E = Dict()
    if kmax == 0
        kmax = maximum(sum(A, dims = 1))
    end
    E = Dict(k => Dict() for k ∈ 2:kmax)

    for j ∈ 1:size(A)[2]
        k = sum(A[:,j])
        if k >= 2
            if k <= kmax
                authors = findall(x -> (x ≈ 1), A[:,j])
                sort!(authors)
                E[k][authors] = get(E[k], authors, 0) + 1
            end
        end
    end

    H = HypergraphModularity.hypergraph(N, E, zero(N))
    HypergraphModularity.computeDegrees!(H);
    return H
end

function readTemporalData(name; kmax = 20, tmin = 0, tmax = 2050)
    
    prefix = "data/$name/$name"
    
    nverts    = open(prefix * "-nverts.txt") do f
        parse.(Int64, collect(eachline(f)))
    end

    simplices = open(prefix * "-simplices.txt") do f
        parse.(Int64, collect(eachline(f)))
    end

    times     = open(prefix * "-times.txt") do f
        parse.(Int64, collect(eachline(f)))
    end
    
    j = 1
    edges = []

    for i ∈ 1:length(nverts)
        e = simplices[j:(j+nverts[i]-1)]
        push!(edges, (e, nverts[i], times[i]))
        j += nverts[i]
    end
    
    subset = [e for e in edges if (length(e[1]) >= 2) & (e[3] ∈ tmin:tmax)];

    nodes = unique(collect(Iterators.flatten([e[1] for e in subset])))
    nodemap = Dict(zip(nodes, 1:length(nodes)));

    function recode(e)
        return [[nodemap[i] for i in e[1]], e[2], e[3]]
    end

    subset = [recode(e) for e in subset];

    
    E = Dict(k => Dict() for k ∈ 2:kmax)

    for e in subset
        sort!(e[1])
        if e[2] ∈ keys(E)
            E[e[2]][e[1]] = get(E[e[2]], e[1], 0) + 1
        end
    end

    H = HypergraphModularity.hypergraph(1:length(nodes), E, [0])
    HypergraphModularity.computeDegrees!(H);
    
    return H
end

