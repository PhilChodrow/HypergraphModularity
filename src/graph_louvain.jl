using Random
using SparseArrays

"""
ConstructAdj: Construct Adjacency List
This takes in a sparse adjacency matrix for a graph, and returns an adjacency
list. While it seems inefficient to store the matrix in multiple ways, as long
as there are no terrible memory issues, it is VERY handy to be able to quickly
get a list of neighbors for each node.

The function also returns the degree vector d. This is NOT the weighted degree,
    d[i] = total number of neighbors of node i
"""
function ConstructAdj(C::SparseMatrixCSC,n::Int64)
    rp = C.rowval
    ci = C.colptr
    Neighbs = Vector{Vector{Int64}}()
    d = zeros(Int64,n)
    for i = 1:n
        # chop up the rp vector and put it in Neighbs
        push!(Neighbs,rp[ci[i]:ci[i+1]-1])
        d[i] = ci[i+1]-ci[i]
    end

    # d is the number of neighbors. This is the unweighted degree,
    # but note importantly that if the original graph is weighted this is
    # not the same as the degree vector d we will sometimes use
    return Neighbs, d
end

"""
Evaluates the LambdaCC objective, which is equivalent to modularity in graphs.

https://dl.acm.org/doi/pdf/10.1145/3178876.3186110

A = adjacency matrix for a graph
c = clustering vector
w = weights vector (often taken to be degrees of nodes)
lam = resolution paramter

This can be replaced with our function for checking the hypergraph modularity
objective.
"""
function LamCCobj(A::SparseMatrixCSC,c,w,lam)

    w_volA = sum(w) # weighted volume
    obj = (lam*(w_volA)^2-lam*sum(w.^2))/2 # constant with respect to the clustering

    for i = 1:maximum(c)
        S = findall(x->x==i,c)
        AS = A[:,S]
        vol = sum(AS.nzval);
        SAS = AS[S,:]
        edges = sum(SAS.nzval);
        cut = vol-edges
        w_volS = sum(w[S])
        obj += 1/2*(cut - lam*w_volS*(w_volA-w_volS))
    end
    return obj
end


"""
The way the Louvain code works, you start with n clusters (singletons), and as nodes
move, many of these become empty. At certain points in the algorithm, it is
helpful to renumber the cluster IDs, so that they go from 1 to M where
M = number of non-empty clusters.

e.g. c = [2 3 5 3 9] ---> c = [1 2 3 2 4]

This version does this with both the Clusters array of arrays AND
the cluster indicator vector c:

c[i] = (integer) cluster ID that node i belongs to
Clusters[j] = (integer array) node IDs for nodes in cluster j
"""
function renumber(c::Vector{Int64},Clusters::Vector{Vector{Int}})

    n = length(c)
    map = unique(c)
    cnew = zeros(Int64,n)

    Clusters = Clusters[map]

    # Rename the clusters
    for i = 1:n
        newClus = findfirst(x->x == c[i],map)
        cnew[i] = newClus
        push!(Clusters[newClus],i)
    end

    return cnew, Clusters

end


"""
See above function, this function doesn't care about the Clusters array
"""
function renumber(c::Vector{Int64})

    n = length(c)
    map = sort(unique(c))
    cnew = zeros(Int64,n)

    # Rename the clusters
    for i = 1:n
        newClus = findfirst(x->x == c[i],map)
        cnew[i] = newClus
    end

    return cnew

end

"""
Run the Louvain algorithm many times, taking the result with the best objective.
"""
function Many_Louvain(A::SparseMatrixCSC{Float64,Int64},w::Vector{Float64},lam::Float64,numtimes::Int64,maxits::Int64=10000)

    n = size(A,1)
    BestObj = Inf
    cBest = collect(1:n)

    for k = 1:numtimes

        ######
        # can be replaced with a "run hypergraph modularity" function
        Cs = LambdaLouvain(A,w,lam,true,maxits)
        ######

        c = Cs[:,end]

        ######
        # can be replaced with a "run hypergraph modularity" function
        obj = LamCCobj(A,c,w,lam)
        ######


        if obj < BestObj
            BestObj = obj
            cBest = c
        end
    end

    return cBest, BestObj
end

"""
The full Louvain algorithm. Lambda is the resolution parameter; think of it as
lambda = gamma/(vol(G)), where gamma is the more traditional resolution parameter
associated with the modularity objective.

w = weights function, w[i] = weight for node i. This is often the degree of node
    i, but doesn't have to be.
"""
function LambdaLouvain(A::SparseMatrixCSC{Float64,Int64},w::Vector{Float64},lam::Float64,randflag::Bool=false,maxits::Int64=10000)

    @assert(issymmetric(A))
    n = size(A,1)

    # Step 1: greedy moves until no more improvement
    #########
    # We can replace this with a greedy moves function for hypergraph modularity
    c, improved = LambdaLouvain_Step(A,w,lam,randflag,maxits)
    #########

    # This code keeps track of the clustering that is found after each call to
    # Step 1.
    if improved
        Cs = c
        c_old = copy(c)
    else
        Cs = c
    end

    # As long as something is still improving each time you call step 1, keep going
    while improved

        # Step 2: Collapse the clustering into supernodes

        ########
        # Should be able to tweak this for hypergraph modularity
        Anew, wnew = collapse_clustering(A,w,c_old)
        # wnew is the weight that comes from summing weights in the same cluster
        # Anew is the reduced adjacency matrix of supernodes
        ########

        # Step 1: Go back to greedy local moves, this time on the reduced graph
        cSuper, improved = LambdaLouvain_Step(Anew,wnew,lam,randflag,maxits)
        N = length(wnew)    # N = number of supernodes = number clusters from last round

        # Undo the procedure that merged nodes into super nodes, so that you get a clustering vector of length n
        # Do this only if the last call to Step 1 led to at least one greedy move.

        if improved

            # Extract what that new clustering is
            c_new = zeros(Int64,n)

            # For each supernode, place all original node IDs that make it up
            # into a cluster of the supernodes label
            for i = 1:N

                # Get the cluster that supernode i is in
                SuperI_Cluster = cSuper[i]

                # Get individual node ID that are in supernode i.

                # findall is slow, but this typically won't need to be called many times
                SuperI_nodes = findall(x->x==i,c_old)
                c_new[SuperI_nodes] .= SuperI_Cluster
            end
            Cs = [Cs c_new]
            c_old = copy(c_new)
        end
    end

    return Cs
end

"""
Run Step 1 of the Louvain algorithm: iterate through nodes and greedily move
nodes to adjacent clusters.
"""
function LambdaLouvain_Step(A::SparseMatrixCSC{Float64,Int64},w::Vector{Float64},lam::Float64,randflag::Bool=false,maxits::Int64=Inf)
    @assert(issymmetric(A))
    n = size(A,1)
    # println("Merging $n Communities")

    # This permutes the node labels, to add randomization in the Louvain
    # algorithm so that you don't always traverse the nodes in the same order
    if randflag
        p = randperm(n)
        A = A[p,p]
        undop = sortperm(p)
        w = w[p]
    end

    # All nodes start in singleton clusters
    c = collect(1:n)
    @assert(size(w,1) == n)
    improving = true
    its = 0

    # This stores the graph as an an adjacency list: it's super helpful
    # to be able to quickly get the set of neighbors for a given node
    # Neighs[i] lists nodes that share an edge with node i
    # degvec[i] = number of neighbors of node i
    Neighbs, degvec = ConstructAdj(A,n)

    # Repeatedly calilng "findall" on the cluster vector c is an inefficient
    # way to get all the nodes that belong in one cluster.
    # Instead, we store a vector of n clusters
    # By the end, many of these will be empty clusters! But it's still
    # faster this way.
    Clusters = Vector{Vector{Int64}}()
    for v = 1:n
        push!(Clusters, Vector{Int64}())
        push!(Clusters[v],v)    # singleton clusters--one for each node
    end

    # As long as we can keep improving, continue
    while improving && its < maxits
        its += 1
        # println("Iteration $its")
        improving = false
        nextclus = 1

        # visit each node in turn
        for i = 1:n

            # println("\n\t NODE $i")

            # Cluster index for node i
            Ci_ind = c[i]

            # Get the indices of nodes in i's cluster
            # Ci = findall(x->x == Ci_ind,c) # findall is slow!
            Ci = Clusters[Ci_ind]

            # Get the indices of i's neighbors--these define clusters
            # we MIGHT move to.
            #   NOT SURE IF THIS WILL BE EXACTLY THE SAME FOR HYPERGRAPH MODULARITY!

            # Ni = findall(x->x>0, A[i,:]) # findall is slow!
            Ni = Neighbs[i]

            ##################################################################
            # THIS IS THE PART THAT REQUIRES THE MOST UPDATING FOR HYPERGRAPH MODULARITY!
            # This is where we compute whether a greedy local move
            # would lead to an improved objective score. The next
            # several blocks of code can hopefully eventually be replaced
            # by a call to a function that does a greedy local move for one node.
            ##################################################################

            ## First few lines here essentially compute the current contribution
            # this node i has on the objective, based on its current cluster placement

                # Weight of negative mistakes currently at i
                neg_inner = w[i]*(sum(w[Ci]) - w[i])

                # Weight of positive mistakes if we remove i from Ci
                pos_inner = sum(A[i,Ci])

                # Increase in mistakes if we remove i from Ci
                total_inner = pos_inner - lam*neg_inner

                BestImprove = 0
                BestC = Ci_ind


            # Get the neighboring clusters of i:
            # there are the only places we would even consider moving node i to
            NC = unique(c[Ni])

            # Now let's see if it's better to move to a nearby cluster, Cj
            for j = 1:length(NC)
                Cj_ind = NC[j]

                # Check how much it would improve to move i to to cluster j
                if Cj_ind == Ci_ind

                    change = 0

                else

                    # Cj = findall(x->x == Cj_ind,c)
                    Cj = Clusters[Cj_ind]

                    # Find the neighbors of i in Cj
                    # cj_neighbs = intersect(Ni,Cj)

                    # Moving i from Ci to Cj adds negative mistakes
                    neg_outer = w[i]*(sum(w[Cj]))

                    # Moving i from Ci to Cj decreases positive mistakes
                    # pos_outer = sum(A[i,cj_neighbs])
                    pos_outer = sum(A[i,Cj])

                    total_outer = lam*neg_outer - pos_outer

                    # This is the overall change in objective if we move i to Cj
                    change = total_outer + total_inner

                end

                # Check if this is currently the best possible greedy move to make
                if change < BestImprove
                    BestImprove = change
                    BestC = Cj_ind
                    improving = true
                end
            end
            # Done checking whether there's a better place to greedily move the cluster to.

            # Move i to the best cluster to move it to, do this
            # only if it strictly improves the objective
            if BestImprove < 0
                ci_old = c[i]
                c[i] = BestC

                # Remove i from its old cluster...
                Clusters[ci_old] = setdiff(Clusters[ci_old],i)

                # ...and add it to its new cluster
                push!(Clusters[BestC],i)

                improving = true # we have a reason to keep iterating!
            end

        end
    end
    if its == 1
        improved = false
    else
        improved = true
        c, Clusters = renumber(c,Clusters)
    end

    if randflag
        c = c[undop]    # if you previously mixed up node order, re-order them now
    end

    return c, improved

end

"""
Step 2 of the Louvain algorithm:
Collapse a clustering into a new network of supernodes and weighted edges, so
that you can then run the same greedy node-moving on supernodes.
"""
function collapse_clustering(A::SparseMatrixCSC{Float64,Int64},w::Vector{Float64},c::Vector{Int64})

    n = size(A,1)
    # Storing clusters in an array of arrays help speed many aspects of the
    # code up.
    # Here we're just setting up singleton clusters for each node.
    Clusters = Vector{Vector{Int64}}()
    for v = 1:maximum(c)
        push!(Clusters, Vector{Int64}())
    end

    for i = 1:n
        push!(Clusters[c[i]],i)
    end

    # Number of supernodes to form = number of clusters
    N = round(Int64,maximum(c))

    I = Vector{Int64}()
    J = Vector{Int64}()
    V = Vector{Float64}()
    wnew = zeros(N)
    Anew = zeros(N,N)
    start = time()

    # Construct a new sparse matrix with new node weights
    for i = 1:N
        #Ci = findall(x->x == i,c)
        # using findall is extremely slow. Here's a faster version of this:
        Ci = Clusters[i]
        wnew[i] = sum(w[Ci])
        ACi = A[Ci,:]
        for j = i+1:N
            Cj = Clusters[j]
            Eij = sum(ACi[:,Cj])
            Anew[i,j] = Eij
        end
    end
    getedges = time()-start

    start = time()
    # Anew = sparse(I,J,V,N,N)
    Anew = Anew+Anew'
    Anew = sparse(Anew)

    return Anew, wnew
end


## Sort sizes: Return a set of cluster sizes, arranged in order
function ClusterSizes(c)
    c = renumber(c)

    sizes = Vector{Int64}()
    for i = 1:maximum(c)
        inds = findall(x->x==i,c)
        push!(sizes,length(inds))
    end

    return sort(sizes,rev=true)
end
