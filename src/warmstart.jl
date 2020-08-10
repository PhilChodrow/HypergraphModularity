# Warm-start
include("graph_louvain.jl")
include("HSBM.jl")

function CliqueExpansion(H::hypergraph,weighted::Bool=true)
    """
    Weighted clique expansion where a hyperedge e is expanded to a
    weighted clique with each edge having weight 1/(|e| - 1)
    """
    n = length(H.D)
    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    ks = setdiff(keys(H.E),1)
    for k in ks
        for edge in keys(H.E[k])
            weight = H.E[k][edge]
            for i = 1:k-1
                ei = edge[i]
                for j = i+1:k
                    ej = edge[j]
                    push!(I,ei)
                    push!(J,ej)
                    if weighted
                        push!(V,weight/(k-1))
                    else
                        push!(V,weight)
                    end
                end
            end
        end
    end
    A = sparse(I,J,V,n,n)
    for i = 1:n; A[i,i] = 0.0; end
    dropzeros!(A)
    A = sparse(A+A')
    return A

end

function CliqueExpansionModularity(H::hypergraph,gamma::Float64=1.0,weighted::Bool=true)
    """
    Perform a clique expansion on the hypergraph H and then run vanilla
    modularity on the resulting graph.
    """
    A = CliqueExpansion(H,weighted)
    return VanillaModularity(A,gamma)
end


function VanillaModularity(A::SparseMatrixCSC{Float64,Int64},gamma::Float64=1.0,randflag::Bool=false,maxits::Int64=10000)
    """
    Vanilla modularity algorithm, obtained by calling the LambdaLouvain algorithm
    implementation from:

    Parameterized Correlation Clustering in Hypergraphs and Bipartite Graphs
    https://arxiv.org/abs/2002.09460

    Code: https://github.com/nveldt/ParamCC/blob/master/src/Graph_Louvain.jl
    """

    d = vec(sum(A,dims = 2))
    vol = sum(d)
    lam = gamma/vol
    Cs = LambdaLouvain(A,d,lam,randflag,maxits)

    c = Cs[:,end]
    @assert(length(c) == n)

    return c
end
