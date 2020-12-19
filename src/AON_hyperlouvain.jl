function SuperNode_PPLouvain(H::hypergraph,Ω::IntensityFunction;α,clusterpenalty=0,kmax=maximum(keys(H.E)),maxits::Int64=100,bigInt::Bool=true,verbose=true,scan_order="random",Z0 = collect(1:length(H.D)))
    randflag = !(scan_order == "lexical")
    cut_weights, vol_weights, e2n, n2e,w,d,elen = AON_Inputs(H,Ω.ω,α,kmax)
    Zset = SuperNode_PPLouvain(n2e,e2n,w,d,elen,cut_weights,vol_weights,kmax,randflag,maxits,verbose,Z0,clusterpenalty);
    Z = Zset[:,end];
    return Z
end

function SuperNode_PPLouvain(H::hypergraph,Ω::IntensityFunction,kmax::Int64 = maximum(keys(H.E)),maxits::Int64=100,bigInt::Bool=true;α,verbose=true,scan_order="random", Z0 = collect(1:length(H.D)))
    Z = SuperNode_PPLouvain(H,Ω;α=α,kmax=kmax,maxits=maxits,bigInt=bigInt,verbose=verbose,scan_order=scan_order,Z0 = Z0)
    return Z
end

function SuperNode_PPLouvain(node2edges::Vector{Vector{Int64}},
    edge2nodes::Vector{Vector{Int64}},w::Vector{Float64},
    d::Vector{Float64},elen::Vector{Int64},
    alp::Vector{Float64},bet::Vector{Float64},
    kmax::Int64,randflag::Bool=false,maxits::Int64=100,verbose::Bool=true,
    Zwarm::Vector{Int64}=Vector{Int64}(),clusterpenalty = 0)

    """
    Supernode version of Hypergraph Louvain for planted partition model
        (all-or-nothing cut penalties).

        * node2edges[i] = list of edge labels that node i is in
        * edge2nodes[e] = list of node labels that edge e contains
        * w[e] = weight of hyperedge e
        * d[i] = degree or weight of node i
        * elen[e] = number of nodes in hyperedge e
                     note that this may be > length(edge2nodes[e])
                     if some of the nodes are supernodes
        * alp = scaling values of for hyperedge cuts
        * bet = scaling values for volumes
        * kmax = maximum hyperedge size
        * Zwarm = warm start clustering (optional)
        * randflag = whether or not to permute node order
        * maxits = Maximum # of greedy passes over the node set
        * clusterpenalty = we include term clusterpenalty*log(k) in objective
            so if clusterpenalty > 0, there is incentive to form fewer clusters
    """

    n = length(d)
    # Step 1: greedy moves until no more improvement
    Z, improved = ANHL_Step(node2edges,edge2nodes,w,d,elen,alp,bet,kmax,
                            randflag,maxits,verbose,Zwarm,clusterpenalty)

    # Store all clusterings found
    if improved
        Zs = Z
        Z_old = copy(Z)
    else
        Zs = Z
    end

    # As long as something is still improving each time you call step 1, keep going
    while improved

        # Step 2: Collapse the clustering into supernodes
        @assert(minimum(Z_old) > 0)
        e2n = deepcopy(edge2nodes)
        Change_id_H!(e2n,Z_old)
        # Merge hyperedges that are the same and remove size 1 hyperedges
        e2n, wSuper, elenSuper =  Condense_H(e2n,w,elen)
        He2n = elist2incidence(e2n,maximum(Z_old))
        n2e = incidence2elist(He2n,true)
        dSuper = condense_d(d,Z_old)
        ########
        # Step 1: Go back to greedy local moves, this time on the reduced hypergraph
        Zsuper, improved = ANHL_Step(n2e,e2n,wSuper,dSuper,elenSuper,alp,bet,
                                kmax,randflag,maxits,verbose,Vector{Int64}(),clusterpenalty)
        N = length(Zsuper)    # N = number of supernodes = number clusters from last round

        @assert(minimum(Zsuper)>0)
        @assert(length(Zsuper)==length(dSuper))
        # println("$N, $improved")
        # Undo the procedure that merged nodes into super nodes, so that you get a clustering vector of length n
        # Do this only if the last call to Step 1 led to at least one greedy move.
        if improved
            # Extract what that new clustering is
            Z_new = zeros(Int64,n)

            # For each supernode, place all original node IDs that make it up
            # into a cluster of the supernodes label
            for i = 1:N

                # Get the cluster that supernode i is in
                SuperI_Cluster = Zsuper[i]

                # Get individual node ID that are in supernode i.
                SuperI_nodes = findall(x->x==i,Z_old)
                Z_new[SuperI_nodes] .= SuperI_Cluster
            end
            @assert(minimum(Z_new)>0)
            Zs = [Zs Z_new]
            Z_old = copy(Z_new)
        end
    end

    return Zs
end


function ANHL_Step(node2edges::Vector{Vector{Int64}},edge2nodes::Vector{Vector{Int64}},
    w::Vector{Float64},d::Vector{Float64},elen::Vector{Int64},
    alp::Vector{Float64},bet::Vector{Float64},kmax::Int64,
    randflag::Bool=false,maxits::Int64=100,verbose::Bool=true,Zwarm::Vector{Int64}=Vector{Int64}(),clusterpenalty=0)
    """
    Basic step Louvain algorithm: iterate through nodes and greedily move
    nodes to adjacent clusters. Does not form supernodes and does not recurse.

    NV: Should work also for degenerate hyperedges, though some bugs are possible

    Features:
    * has weight for size of hyperedge rather than computing the weights
    * allows you to warm start with a clustering
    * allows you to randomize node order
    """
    if verbose println("One step of all-or-nothing HyperLouvain") end
    m = length(edge2nodes)
    n = length(node2edges)

    # Randomize node visit order
    if randflag
        node_order = Random.randperm(n)
        # NOTE: may be faster to re-arrange
        #   all the data structures so that we march through
        #   Z in order.
    else
        node_order = 1:n
    end

    if verbose println("") end

    # Warm start clustering
    if length(Zwarm) > 0
        Z = renumber(Zwarm)
    else
        # Default is to put all nodes in their own cluster
        Z = collect(1:n)
    end


	K = maximum(Z)   			# number of clusters
    ClusVol = zeros(K)			# volume of clusters
	ClusSize = zeros(K)			# size of each cluster
	for i = 1:n
		ClusVol[Z[i]] += d[i]
		ClusSize[Z[i]] += 1
	end

    # Store cut penalties for each edge
    cutpenalty = zeros(m)
    for e = 1:m
        k = elen[e]      				# size of the edge
        cutpenalty[e] = alp[k]*w[e] 	# penalty for cutting it
    end

    # Store node neighbors of each node
    Neighbs = NeighborList(node2edges, edge2nodes)

    # Store all the edges for the node, minus that node id
    # This makes computing cut changes faster
    edgelists = Vector{Vector{Vector{Int64}}}()
    for i = 1:n
        iedges = Vector{Vector{Int64}}()    # each node has a list of hyperedges it's in, i.e. a list of lists of nodes.
        for e in node2edges[i]
            edge = edge2nodes[e]
            push!(iedges, setdiff(edge,i))  # don't save i itself
        end
        push!(edgelists,iedges)
    end

    improving = true
    iter = 0
    changemade = false

    # Initialize vols and cuts to track efficient local changes
    r = kmax
    toler = 1e-8

    mainstart = time()
    while improving && iter < maxits

        iter += 1
        if mod(iter,1) == 0
            if verbose println("Louvain Iteration $iter") end
        end
        improving = false

        # visit each node in turn
        for i = node_order

            # Cluster index for node i
            Ci_ind = Z[i]

            # Get the size of the cluster
			clus_size = ClusSize[Ci_ind]

            # Get the indices of i's neighbors--these define clusters we might move to
            Ni = Neighbs[i]

            # Get the neighboring clusters of i
            NC = unique(Z[Ni])

            # The default is to not move the node
            BestC = Z[i]
            BestImprove = 0

            # Set of edges that i is in
            Cv = node2edges[i]

            # Set of hyperedges that node i is in, not including i itself
            Cv_list = edgelists[i]

            # Volume of the set currently
            vS = ClusVol[Ci_ind]
            dv = d[i]

            for j = 1:length(NC)

				# Check how much it would improve to move i to to cluster j
                Cj_ind = NC[j]

                if Cj_ind == Ci_ind
                    continue
                end

                # Change in volume
                vJ = ClusVol[Cj_ind]
                Δvol = 0
                for k = 1:kmax
					# Better if this is smaller
                    Δvol += bet[k]*((vS-dv)^k + (vJ+dv)^k - vS^k - vJ^k)
                end

                # Change in cut
                Δcut = 0
                for eid = 1:length(Cv)
                    # change in cut for edge e when v moves from
                    # cluster Ci to cluster Cj
                    e = Cv[eid]
                    edge_noi = Cv_list[eid]

                    if elen[e] > 1
						we = cutpenalty[e]
                        Δcut += move_cut(i,Z,edge_noi,Ci_ind,Cj_ind,we)
                    end
                end

                # Change in objective due to change in cluster number
                Δclus = 0
                if clusterpenalty > 0 && clus_size == 1
                    Δclus = clusterpenalty*(log(K-1)-log(K))
                end
                change = Δcut + Δvol + Δclus # want this to be negative

                # Check if this is currently the best possible greedy move to make
                # The tolerance helps with making some behavior more stable:
                #   you only move if there's a numerically nonzero reason to move.
                if change < BestImprove - toler
                    BestImprove = change
                    BestC = Cj_ind
                    improving = true
                end
            end

            # Move i to the best new cluster, only if it strictly improves modularity
            if BestImprove < -toler

                # update clustering
                ci_old = Z[i]
                Z[i] = BestC

                # Remove i from its old cluster...
                ClusVol[ci_old] -= dv
				ClusSize[ci_old] -= 1

                # ...and add it to its new cluster
                ClusVol[BestC] += dv
				ClusSize[BestC] += 1
                changemade = true
                improving = true # we have a reason to keep iterating!

                if clus_size == 1
                    # we moved a singleton to another cluster
                    K -= 1
                end

            end
        end
    end
    mainloop = time()-mainstart
    if ~changemade
        improved = false
        if verbose println("No nodes moved clusters") end
    else
        improved = true
    end
    if verbose println("Main loop took $mainloop seconds") end
	Z = renumber(Z)
    return Z, improved

end

"""
Compute the change from moving node v from Ci to Cj
"""
function move_cut(v::Int64,Z::Vector{Int64},erest::Vector{Int64},Ci_ind::Int64,
    Cj_ind::Int64,we::Float64)

    # p = Z[erest]            # set of clusters in e, not counting e
    p1 = Z[erest[1]]
    na = notsame(erest,Z,p1)


    if na
        # if the edge is cut regardless of whether v is moved, return 0
        return 0
    else
        # p1 = p[1]   # this the the cluster all other nodes are in
        if Ci_ind == p1
            # then moving v cuts the edge
            return we
        elseif Cj_ind == p1
            # then moving v to Cj_ind "uncuts" the edge
            return -we
        else
            # moving doesn't matter
            return 0
        end
    end
end

"""
e is the set of nodes in the hyperedge. Z stores their clusters
"""
function notsame(e::Vector{Int64},Z::Vector{Int64},p1::Int64)
    for i = 2:length(e)
        if Z[e[i]] != p1
            return true
        end
    end
    return false
end


function Condense_H(edge2nodes::Vector{Vector{Int64}},w::Vector{Float64})
    """
    Condense Hypergraph.

        A hypergraph stored as an edgelist may contain some edges
        that are identical. This function finds them and condenses
        them instead into a weighted hypergraph edgelist
    """
    Hdict = Dict{Vector{Int64},Float64}()   # maps from set of nodes to a weight

    for e = 1:length(edge2nodes)
        edge = edge2nodes[e]
        sort!(edge)
        we = w[e]
        wt = get(Hdict,edge,0)
        Hdict[edge] = wt+we
    end

    newlist = Vector{Vector{Int64}}()
    for edge in keys(Hdict)
        push!(newlist,edge)
        push!(wnew,Hdict[edge])
    end

    return newlist, wnew
end

function Condense_H(edge2nodes::Vector{Vector{Int64}},w::Vector{Float64},label::Vector{Int64})
    """
    Condense Hypergraph.

        A hypergraph stored as an edgelist may contain some edges
        that are identical. This function finds them and condenses
        them instead into a weighted hypergraph edgelist.

        In this version, each hyperedge is additionally associated with a
        label, and we only merge nodes if they have the same node set
        AND the same label.

        The edge label is mainly meant to store original hyperedge size.
        Since some of the nodes may be supernodes, we need to keep track of
        the original hyperedge size.

        Also, get rid of size 1 hyperedges in the process
    """
    Hdict = Dict{Vector{Int64},Float64}()   # maps from set of nodes to a weight

    @assert(length(label)==length(edge2nodes))

    for e = 1:length(edge2nodes)
        edge = edge2nodes[e]
        lab = label[e]
        sort!(edge)
        ky = [lab; edge]    # the key is the label + edge
        we = w[e]
        wt = get(Hdict,ky,0)
        Hdict[ky] = wt+we
    end

    newlist = Vector{Vector{Int64}}()
    wnew = Vector{Float64}()
    labelsnew = Vector{Int64}() # indices of edges you keep
    for ky in keys(Hdict)
        edge = ky[2:end]
        lab = ky[1]
        if length(edge) > 1
            push!(newlist,edge)
            push!(wnew,Hdict[ky])
            push!(labelsnew,lab)
        end
    end

    return newlist, wnew, labelsnew
end


function Change_id_H!(edge2nodes::Vector{Vector{Int64}},Z::Vector{Int64},uniqueflag::Bool=true)
"""
Given a clustering vector Z, change node IDs in the hypergraph
so that node i becomes node Z[i]. This will potentially lead
to many degenerate hyperedges where the same node shows up multiple times.

This is useful as a super-node step in clustering heuristics like louvain.

If uniqueflag == true, don't store degenerate hyperedges. Instead, just save
one copy of the vector
"""

    for e = 1:length(edge2nodes)
        for i = 1:length(edge2nodes[e])
            edge2nodes[e][i] = Z[edge2nodes[e][i]]
        end
        if uniqueflag
            unique!(edge2nodes[e])
        end
    end

    return edge2nodes
end

function condense_d(d::Vector{Float64},Z::Vector{Int64})
    """
    Condense a degree vector to a weights vector by summing up volumes of
    clusters.
    """

    n = maximum(Z)
    dnew = zeros(n)
    for i = 1:n
        S = findall(x->x==i,Z)
        dnew[i] = sum(d[S])
    end
    return dnew
end


function AON_Inputs(H,ω,α,kmax)
    """
    Given a hyperedge H, ω function, and α which parameterizes ω,
    return all of the data structures you need to run the new
    all-or-nothing louvain algorithm
    """
    cut_weights = zeros(kmax)
    vol_weights = zeros(kmax)
    for k = 2:kmax
        cut_weights[k] = log(ω([1,k],α))-log(ω([0,k],α))
        vol_weights[k] = ω([1,k],α)-ω([0,k],α)
    end

    He2n, edge_weights = hypergraph2incidence(H)
    e2n = incidence2elist(He2n);
    n2e = incidence2elist(SparseArrays.sparse(He2n'));
    m = length(e2n)
    edge_len = zeros(Int64,m)
    for e = 1:m
        edge_len[e] = length(e2n[e])
    end
    deg = vec(sum(He2n,dims = 1))
    # @show H.D-deg
    # @assert(deg == H.D)

    return cut_weights, vol_weights, e2n, n2e, edge_weights,deg,edge_len

end



function Elist_to_Hypergraph(elist::Vector{Vector{Int64}}, maxsize::Int64=25)
    """
    Converts a binary hypergraph incidence matrix to type "hypergraph".
    Assumes edges are properly sorted already
    """
    E = Dict{Integer, Dict}()
    for edge in elist
        # sort!(edge)
        if length(edge) > maxsize; continue; end
        sz = length(edge)
        if !haskey(E, sz)
            E[sz] = Dict{}()
        end
        E[sz][edge] = 1
    end

    D = zeros(Int64, n)
    for (sz, edges) in E
        for (e, _) in edges
            D[e] .+= 1
        end
    end
    N = 1:n
    return hypergraph(N, E, D)
end

function Hmat_to_Hypergraph(H::SparseArrays.SparseMatrixCSC, maxsize::Int64=25)
    """
    Converts a binary hypergraph incidence matrix to type "hypergraph"
    """
    E = Dict{Integer, Dict}()
    elist = incidence2elist(H)
    for edge in elist
        sort!(edge)
        if length(edge) > maxsize; continue; end
        sz = length(edge)
        if !haskey(E, sz)
            E[sz] = Dict{}()
        end
        E[sz][edge] = 1
    end

    D = zeros(Int64, n)
    for (sz, edges) in E
        for (e, _) in edges
            D[e] .+= 1
        end
    end
    N = 1:n
    return hypergraph(N, E, D)
end
