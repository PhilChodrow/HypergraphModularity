# simple functions for assisting data analysis

function downSampleEdges!(H::hypergraph, prop)
    """
    prop: the proportion of edges in H to maintain
    """
    for k in keys(H.E)
        for e in keys(H.E[k])
            for j in 1:H.E[k][e]
                if rand() > prop
                    H.E[k][e] -= 1 
                end
            end
            if H.E[k][e] == 0
                delete!(H.E[k], e)
            end
        end
    end
    HypergraphModularity.computeDegrees!(H);
end

function removeDegreeZeroNodes!(H::hypergraph)
    to_remove = H.D .== 0
    println("not implemented")
end

# Create a subhypergraph
# Returns (subhypergraph, node_map) where
#     node_map[i] is the new index of node i
#     the keys of node_map are the indices of in_subhypergraph that evaluate to true
function subhypergraph(h::hypergraph, in_subhypergraph::Vector{Bool})
    # Get new set of edges
    new_edges = []
    for (sz, edges) in h.E
        for (edge, val) in edges
            new_edge = filter(v -> in_subhypergraph[v], edge)
            if length(new_edge) > 1
               push!(new_edges, new_edge)
            end
        end
    end
    
    # renumbering
    node_map = Dict{Int64,Int64}()
    for (i, val) in enumerate(in_subhypergraph)
        if val
            node_map[i] = length(node_map) + 1
        end
    end

    renumber_edge(e) = [node_map[v] for v in e]
    renumbered_new_edges = [renumber_edge(e) for e in new_edges]

    # New edges
    subE = Dict{Integer, Dict}()
    for edge in renumbered_new_edges
        sz = length(edge)
        if !haskey(subE, sz)
            subE[sz] = Dict{}()
        end
        subE[sz][edge] = 1
    end

    # New degrees
    n = length(node_map)
    subD = zeros(Int64, n)
    for (sz, edges) in subE
        for (e, _) in edges
            subD[e] .+= 1
        end
    end
    
    return hypergraph(1:n, subE, subD), node_map
end

function mutualInformation(Z, Ẑ, normalized = false)
    """
    Mutual information between two clusterings, optionally normalized
    Probably can be computed MUCH faster, but unlikely to be a bottleneck 
    in context. 
    """
    n = length(Z)
    
    p_XY = Dict()
    p_X  = Dict()
    p_Y  = Dict()

    for i = 1:length(Z)
        p_XY[(Z[i], Ẑ[i])] = get!(p_XY, (Z[i], Ẑ[i]), 0) + 1/n
        p_X[Z[i]]          = get!(p_X, Z[i], 0)          + 1/n
        p_Y[Ẑ[i]]          = get!(p_Y, Ẑ[i], 0)          + 1/n
    end
    
    I = 0
    for x in keys(p_X), y in keys(p_Y)
        try
            I += p_XY[x,y]*log(p_XY[x,y]/(p_X[x]*p_Y[y]))
        catch e
            nothing
        end
    end
    
    if normalized
        H_X, H_Y = 0, 0
        for x in keys(p_X)
            H_X -= log(p_X[x])*p_X[x]
        end
        for y in keys(p_Y)
            H_Y -= log(p_Y[y])*p_Y[y]
        end
        return (2*I)/(H_X + H_Y)
    end
    
    return I
end


function subHypergraph(H, b, Z = nothing)
    """
    b a boolean vector of the same length as H.N, indicating
    which nodes should be included
    VERY SLOW at the moment, could likely be improved. 
    """
    key = findall(x -> x == 1, b)
    nodemap = Dict(zip(key, 1:length(key)))
    
    N_ = copy(H.N)
    
    ix = [N_[i] ∈ key for i in 1:length(N_)]
    N_ = N_[ix]
    
    if !isnothing(Z)
        Z_ = copy(Z)
        Z_ = Z_[ix]
    end
    
    N_ = [nodemap[i] for i in N_]
    
    E_ = copy(H.E)
    for k in keys(E_)
        filter!(e -> all(j in key for j in first(e)), E_[k])
    end
    
    E__ = Dict(k => Dict() for k ∈ 2:maximum(keys(H.E)))
    for k in keys(E__)
        for e in keys(E_[k])
            e_ = [nodemap[i] for i in e]
            sort!(e_)
            E__[k][e_] = get(E__[k],e_, 0) + 1
        end
    end
    
    H_ = hypergraph(N_, E__, [])
    HypergraphModularity.computeDegrees!(H_)
    if !isnothing(Z)
        return H_, Z_
    end
    return H_
end

function projectedGraph(H)
    
    n = length(H.D)
    A = CliqueExpansion(H, false, false)
    ix, jx, w = SparseArrays.findnz(A)
    E = Dict(sort([ix[k], jx[k]]) => w[k] for k in 1:length(ix))

    H̄ = hypergraph(collect(1:n), Dict(1 => Dict(), 2 => E), [0])
    computeDegrees!(H̄);
    return H̄
end

function kcore(H, Z, core)
    H_, Z_ = copy(H), copy(Z)
    for i in 1:10
        H_, Z_ = subHypergraph(H_, H_.D .>= core, Z_)
    end
    return H_, Z_
end

function removeEdges!(H; remove = [])
    for k in remove
        pop!(H.E, k, nothing)
    end
    computeDegrees!(H)
end