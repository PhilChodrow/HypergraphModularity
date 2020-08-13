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
