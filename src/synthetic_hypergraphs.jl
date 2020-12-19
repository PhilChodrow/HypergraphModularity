using StatsBase
using Clustering
using LinearAlgebra


function GenerateHypergraphAll(n,m,K,pvals,rmin,rmax,cluster_sizes=ones(K),r_sizes=ones(rmax-rmin),cluster_prefs=ones(K))
    """
    Generate a synthetic hypergraph, store in incidence matrix
    Output several other things about hypergraph like edge sizes and degrees
    and edge list.
    * Set a fixed number of hyperedges m and clusters K (you can randomize this outside this function)
    * Set a distribution determining the relative proportion of hyperedges of each size from r_min to r_max
    * For each hyperedge, with probability p, assign it to one community at random,
        selecting a r nodes uniformly at random for a hyperedge of size r
    * Otherwise, with probability 1-p, select r nodes uniformly at random from across the entire hypergraph
    """

    # Normalize to get proportions/probabilities
    r_sizes = Weights(r_sizes/sum(r_sizes))
    cluster_sizes = Weights(cluster_sizes/sum(cluster_sizes))         # proportion of nodes in each cluster
    cluster_prefs = Weights(cluster_prefs/sum(cluster_prefs))        # relative proportion of interior hyperedge in each cluster

    # Generate the ground truth clustering of nodes
    # ground_truth = GenerateTruth(n,cluster_sizes)
    Clusters = Vector{Vector{Int64}}()
    for j = 1:K
        push!(Clusters,Vector{Int64}())
        @assert(n*cluster_sizes[j] > 2*rmax)
    end


    ground_truth = zeros(Int64,n)
    for i = 1:n
        c = sample(1:K,cluster_sizes)
        ground_truth[i] = c
        push!(Clusters[c],i)
    end


    U = Vector{Int64}()
    E = Vector{Int64}()
    EdgeList = Vector{Vector{Int64}}()
    E_lengths = zeros(Int64,m)
    for enum = 1:m
        e = SampleEdge(pvals,n,rmin,rmax,cluster_prefs,r_sizes,Clusters)
        push!(EdgeList,e)
        E_lengths[enum] = length(e)
        for node in e
            push!(U,node)
            push!(E,enum)
        end
    end
    He2n = SparseArrays.sparse(E,U,ones(length(U)),m,n)
    deg = vec(sum(He2n,dims = 1))
    return He2n, EdgeList, E_lengths, deg, ground_truth

end

function GenerateHypergraph(n,m,K,p,rmin,rmax,cluster_sizes=ones(K),r_sizes=ones(rmax-rmin),cluster_prefs=ones(K))
    """
    Generate a synthetic hypergraph in the following way
    * Set a fixed number of hyperedges m and clusters K (you can randomize this outside this function)
    * Set a distribution determining the relative proportion of hyperedges of each size from r_min to r_max
    * For each hyperedge, with probability p, assign it to one community at random,
        selecting a r nodes uniformly at random for a hyperedge of size r
    * Otherwise, with probability 1-p, select r nodes uniformly at random from across the entire hypergraph
    """

    # Normalize to get proportions/probabilities
    r_sizes = Weights(r_sizes/sum(r_sizes))
    cluster_sizes = Weights(cluster_sizes/sum(cluster_sizes))         # proportion of nodes in each cluster
    cluster_prefs = Weights(cluster_prefs/sum(cluster_prefs))        # relative proportion of interior hyperedge in each cluster

    # Generate the ground truth clustering of nodes
    # ground_truth = GenerateTruth(n,cluster_sizes)
    Clusters = Vector{Vector{Int64}}()
    for j = 1:K
        push!(Clusters,Vector{Int64}())
        @assert(n*cluster_sizes[j] > 2*rmax)
    end


    ground_truth = zeros(Int64,n)
    for i = 1:n
        c = sample(1:K,cluster_sizes)
        ground_truth[i] = c
        push!(Clusters[c],i)
    end

    EdgeList = Vector{Vector{Int64}}()
    for i = 1:m
        edge = SampleEdge(p,n,rmin,rmax,cluster_prefs,r_sizes,Clusters)
        push!(EdgeList,edge)
    end

    return EdgeList, ground_truth

end

function SampleEdge(pvals,n,rmin,rmax,cluster_prefs,r_sizes,Clusters)
    replace = false
    r = sample(rmin:rmax,r_sizes)       # select a hyperedge size
    if rand(1)[1] < pvals[r-rmin+1]

        c = sample(1:K,cluster_prefs)   # select a random cluster to put the edge in
        return sort(sample(Clusters[c],r,replace = replace))
    else
        return sort(sample(1:n,r,replace = replace))
    end
end


function ari(x,y)
    evaluations = randindex(x, y)
    ari = evaluations[1]
    return ari
end



function nmi(x::Vector{Int64}, y::Vector{Int64})
# Compute normalized mutual information I(x,y)/sqrt(H(x)*H(y)) of two discrete variables x and y.
# Input:
#   x, y: two integer vector of the same length
# Ouput:
#   z: normalized mutual information z=I(x,y)/sqrt(H(x)*H(y))


@assert(length(x) == length(y));
n = length(x)
# x = reshape(x,1,n);
# y = reshape(y,1,n);

l = minimum([minimum(x),minimum(y)])
x = x.-(l-1)
y = y.-(l-1)
k = maximum([maximum(x),maximum(y)])

idx = collect(1:n)
Mx = sparse(idx,x,ones(length(x)),n,k,n)
My = sparse(idx,y,ones(length(x)),n,k,n)
MxMyn = Mx'*My/n
Pxy = MxMyn.nzval    #joint distribution of x and y
Hxy = -dot(Pxy,log2.(Pxy))


# hacking, to elimative the 0log0 issue
Px = Vector{Float64}()
Py = Vector{Float64}()

MxMean = mean(Mx,dims=1)
MyMean = mean(My,dims=1)
for i = 1:length(MxMean)
    if MxMean[i] > 0
        push!(Px,MxMean[i])
    end
end
for i = 1:length(MyMean)
    if MyMean[i] > 0
        push!(Py,MyMean[i])
    end
end

# entropy of Py and Px
Hx = -dot(Px,log2.(Px));
Hy = -dot(Py,log2.(Py));

# mutual information
MI = Hx + Hy - Hxy

# normalized mutual information
z = sqrt((MI/Hx)*(MI/Hy))
z = maximum([0,z])

return z
end


function learn_alpha_wrapper(H,c,kmax,n)
    """
    Given a hypergraph and clustering,
    learn the alpha parameters and put them in the right format
    for AON_Louvain.jl
    """
    α = zeros(2*kmax)
    function ω(p,α)
        k = p[2]
        δ = p[1]
        return ((1+(1-δ))*n)^α[k] / (n^α[k + kmax])
    end
    Ω = allOrNothingIntensityFunction(ω, kmax)
    st = time()
    α = learnParameters(H, c, Ω, α; n_iters = 10, amin = -10, amax = 10)
    learntime = time()-st
    cut_weights = zeros(kmax)
    vol_weights = zeros(kmax)
    for k = 2:kmax
        cut_weights[k] = log(ω([1,k],α))-log(ω([0,k],α))
        vol_weights[k] = ω([1,k],α)-ω([0,k],α)
    end
    cut_weights, vol_weights, learntime, α
end
