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
