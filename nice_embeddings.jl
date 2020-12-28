using UMAP

using Pkg
Pkg.activate(".")
using HypergraphModularity

using Arpack
using BSON
using LinearAlgebra
using SparseArrays
using Random

function regularized_spectral_embedding(A; tol=1e-10, maxiter=300, nev=3)
    dC = vec(sum(A,dims=1))
    dR = vec(sum(A,dims=2))
    τC = sum(dC) / length(dC)
    τR = sum(dR) / length(dR)
    DC = Diagonal(1.0 ./ sqrt.(dC .+ τC))
    DR = Diagonal(1.0 ./ sqrt.(dR .+ τR))
    N = DR * A * DC
    s = svds(N, nsv=nev, tol=tol, maxiter=maxiter)
    U, V = s[1].U, s[1].V
    return U, V
end

function read_dataset(dataset)
    H, labels = read_hypergraph_data(dataset, 25, 2)
    names = read_hypergraph_label_names(dataset)
    II = Int64[]
    JJ = Int64[]
    edge_ind = 1
    for E in values(H.E), (edge, _) in E
        for v in edge
            push!(II, v)
            push!(JJ, edge_ind)
	end
        edge_ind += 1
    end
    B = sparse(II, JJ, 1, maximum(II), maximum(JJ));
    size(B)
    return H, labels, names, B
end

function gen_embeddings(dataset)
    H, labels, names, B = read_dataset(dataset)
    U, V = regularized_spectral_embedding(B, tol=1e-10, maxiter=500, nev=25)
    xy_coords = umap([U' V'], 2)
    bson("$dataset-xy.bson", Dict("xy" => xy_coords, "U" => U, "V" => V, "B" => B, "labels" => labels, "names" => names))
end
