# Warm-start
include("graph_louvain.jl")
include("HSBM.jl")

function CliqueExpansion(H::hypergraph,weighted::Bool=true,binary::Bool=false)
    """
    Weighted clique expansion where a hyperedge e is expanded to a
    weighted clique with each edge having weight 1/(|e| - 1)
    """
    n = length(H.D)
    I = Int64[]
    J = Int64[]
    V = Float64[]
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
                        push!(V, weight / (k - 1))
                    else
                        push!(V, weight)
                    end
                end
            end
        end
    end
    A = SparseArrays.sparse(I,J,V,n,n)
    for i = 1:n
        A[i, i] = 0.0
    end
    SparseArrays.dropzeros!(A)
    A = SparseArrays.sparse(A+A')
    if binary
        I, J, V = SparseArrays.findnz(A)
        A = SparseArrays.sparse(I, J, 1, n, n)
    end
    return A
end

function CliqueExpansionModularity(H::hypergraph,gamma::Float64=1.0,weighted::Bool=true,randflag::Bool=false,binary::Bool=false)
    """
    Perform a clique expansion on the hypergraph H and then run vanilla
    modularity on the resulting graph.
    """
    A = CliqueExpansion(H,weighted,binary)
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



function computeDyadicResolutionParameter(H, Z; mode = "γ", weighted=true, binary=false)
    """
    compute the dyadic resolution parameter associated to a partition using the formula from Newman (2016): https://arxiv.org/abs/1606.02319
    """

    G = CliqueExpansion(H, weighted, binary)
    I, J = SparseArrays.findnz(G)
    n = maximum(I)  # number of nodes
    m = sum(G)/2    # number of edges

    # form degree sequence and edge counts
    D = vec(sum(G, dims=1))

    m_in = 0
    m_out = 0

    for k in 1:length(I)
        if Z[I[k]] == Z[J[k]]
            m_in  += G[I[k], J[k]]/2
        else
            m_out += G[I[k], J[k]]/2
        end
    end

    # compute resolution parameter
    V = [sum(D[Z .== c]) for c in unique(Z)]
    ω_in = 4*m*m_in / (sum(V.^2))
    ω_out = (2m - 2m_in)/(2m - (sum(V.^2)/(2m)))
    
    if mode == "γ"
        γ = (ω_in - ω_out)/(log(ω_in) - log(ω_out))
        return γ
    else
        return(ω_in, ω_out)
    end
end


function dyadicModularity(H, Z, γ; weighted=true, binary=false)
    G = CliqueExpansion(H, weighted, binary)
    d = vec(sum(G, dims=1))

    # non-degree (cut) term
    edge_obj = 0.0
    for (i, j, v) in zip(SparseArrays.findnz(G)...)
        if Z[i] == Z[j]
            edge_obj += v
        end
    end

    # volume terms
    vols = Dict{Int64, Float64}()
    for c in unique(Z)
        vols[c] = 0.0
    end
    for i = 1:length(d)
        vols[Z[i]] += d[i]
    end

    Q = edge_obj
    volG = sum(d)
    vol_term = 0.0
    for c in unique(Z)
        Q -= γ * vols[c]^2 / volG
    end

    return Q / volG
end

function dyadicLogLikelihood(H, Z, ω_in, ω_out; weighted=false, binary=false)
    G = CliqueExpansion(H, weighted, binary)
    d = vec(sum(G, dims=1))

    # Eq. (14) from https://arxiv.org/pdf/1606.02319.pdf
    γ = (ω_in - ω_out) / (log(ω_in) - log(ω_out))
    Q = dyadicModularity(H, Z, γ; weighted=weighted)
    m = sum(d) / 2
    B = m * log(ω_in / ω_out)
    C = m * (ω_out + log(ω_out))
    return  B * Q + C
end
