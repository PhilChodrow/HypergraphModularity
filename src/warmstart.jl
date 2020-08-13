# Warm-start
include("graph_louvain.jl")
include("HSBM.jl")

function CliqueExpansion(H::hypergraph,weighted::Bool=true,binary::Bool=false)
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
    A = SparseArrays.sparse(I,J,V,n,n)
    for i = 1:n; A[i,i] = 0.0; end
    SparseArrays.dropzeros!(A)
    A = SparseArrays.sparse(A+A')
    if binary
        I,J,V = findnz(A)
        A = sparse(I,J,ones(length(I)),n,n)
    end
    return A
end

function CliqueExpansionModularity(H::hypergraph,gamma::Float64=1.0,weighted::Bool=true,randflag::Bool=false)
    """
    Perform a clique expansion on the hypergraph H and then run vanilla
    modularity on the resulting graph.
    """
    A = CliqueExpansion(H,weighted)
    return VanillaModularity(A,gamma,randflag)
end


function VanillaModularity(A::SparseArrays.SparseMatrixCSC{Float64,Int64},gamma::Float64=1.0,randflag::Bool=false,maxits::Int64=10000)
    """
    Vanilla modularity algorithm, obtained by calling the LambdaLouvain algorithm
    implementation from:

    Parameterized Correlation Clustering in Hypergraphs and Bipartite Graphs
    https://arxiv.org/abs/2002.09460

    Code: https://github.com/nveldt/ParamCC/blob/master/src/Graph_Louvain.jl
    """


    d = vec(sum(A,dims = 2))
    n = length(d)
    vol = sum(d)
    lam = gamma/vol
    Cs = LambdaLouvain(A,d,lam,randflag,maxits)

    c = Cs[:,end]
    @assert(length(c) == n)

    return c
end



function computeDyadicResolutionParameter(H, Z)
    """
    compute the dyadic resolution parameter associated to a partition using the formula from Newman (2016): https://arxiv.org/abs/1606.02319
    """

    G = CliqueExpansion(H)
    I, J = SparseArrays.findnz(G)
    n = maximum(I) # number of nodes
    m = sum(G)/2     # number of edges

    # form degree sequence and edge counts
    D = zero(1.0.*collect(1:n))

    m_in = 0
    m_out = 0

    for k in 1:length(I)
        D[I[k]] += G[I[k], J[k]]
        if Z[I[k]] == Z[J[k]]
            m_in  += G[I[k], J[k]]/2
        else
            m_out += G[I[k], J[k]]/2
        end
    end

    # compute resolution parameter
    V = [sum(D[Z.==k]) for k in 1:maximum(Z)]
    ω_in = 4*m*m_in / (sum(V.^2))
    ω_out = (2m - 2m_in)/(2m - (sum(V.^2)/(2m)))
    γ = (ω_in - ω_out)/(log(ω_in) - log(ω_out))

    return(γ)
end

function dyadicModularity(H, Z, γ)
    G = CliqueExpansion(H)
    I, J = SparseArrays.findnz(G)
    n = maximum(I) # number of nodes
    m = sum(G)/2     # number of edges

    D = zero(1.0.*collect(1:n))

    for k in 1:length(I)
        D[I[k]] += G[I[k], J[k]]
    end

    Q = 0.0

    for i = 1:n, j = 1:n
        if Z[i] == Z[j]
            Q += G[i, j] - γ*(D[i]*D[j])/(2m)
        end
    end
    Q = Q/(2m)
end
