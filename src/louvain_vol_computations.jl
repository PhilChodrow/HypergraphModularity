include("omega.jl")
include("HSBM.jl")
include("cut.jl")
include("vol.jl")
include("utils.jl")
include("dict_ops.jl")
include("diffs.jl")
include("inference.jl")
include("objectives.jl")
include("graph_louvain.jl")
include("warmstart.jl")
include("hyper_format.jl")
include("hyperlouvain_helpers.jl")
include("read_data.jl")

function HyperLouvain_Vols(H::hypergraph,kmax::Int64,alp,bet,Ω::IntensityFunction,checkvols::Bool=false,maxits::Int64=100,bigInt::Bool=true;α,verbose=true,scan_order="lexical")
    """
    Basic step Louvain algorithm: iterate through nodes and greedily move
    nodes to adjacent clusters. Does not form supernodes and does not recurse.

    H: hypergraph
    Ω: group interaction function, as constructed by ΩFromDict(D)
    bigInt::Bool: whether to convert the degree sequence to an array of BigInt when evaluating the volume term in second_term_eval(). Recommended.
    return: Z::array{Int64, 1}: array of group labels.
    """
    He2n, w = hypergraph2incidence(H)
    edge2nodes = incidence2elist(He2n);
    Hyp, w = hyperedge_formatting(H)    # hyperedge to node list
    node2edges = EdgeMap(H)             # node to hyperedge list
    n = length(H.D)

    Neighbs = NeighborList(node2edges, Hyp)  # Store neighbors of each node


    d = vec(sum(He2n,dims = 1))

    kmin = minimum(keys(H.E))
    Z = collect(1:n)                    # All nodes start in singleton clusters
    Clusters = Vector{Vector{Int64}}()  # Also store as cluster list
    for v = 1:n
        push!(Clusters, Vector{Int64}())
        push!(Clusters[v],v)
    end

    iter = 0
    improving = true
    changemade = false

    # Data structures for fast local changes in volumes
    r = kmax
    V, μ, M = evalSums(Z, H.D, r;constants=false, bigInt=bigInt)
    C = evalConstants(r)
    cuts = evalCuts(H, Z, Ω)

    # initialize memoization of edge penalties
    Pn           = [log(Ω.ω(Ω.P(collect(1:i)), α)) for i in 1:kmax]
    edge2penalty = [w[e]*Pn[length(Hyp[e])] for e in 1:length(Hyp)]

    m = length(Hyp)
    elen = zeros(Int64,m)
    for e = 1:m
        elen[e] = length(edge2nodes[e])
    end

    # Save hyperedges in a special way for easier cut checking:
    edgelists = Vector{Vector{Vector{Int64}}}()
    for i = 1:n
        iedges = Vector{Vector{Int64}}()    # each node has a list of hyperedges it's in, i.e. a list of lists of nodes.
        for e in node2edges[i]
            edge = edge2nodes[e]
            push!(iedges, setdiff(edge,i))  # don't save i itself
        end
        push!(edgelists,iedges)
    end


    tic = time()
    # Main Greedy Local Move Loop
    while improving && iter < maxits

        iter += 1
        if mod(iter,1) == 0;
            if verbose println("Louvain Iteration $iter") end
        end
        improving = false


        scan = scan_order == "lexical" ? (1:n) : Random.shuffle(1:n)

        for i = scan                     # visit each node in turn

            Ci_ind = Z[i]               # Cluster index for node i
            Ci = Clusters[Ci_ind]       # Indices of nodes in i's cluster
            Ni = Neighbs[i]             # Get the indices of i's neighbors
            NC = unique(Z[Ni])          # Clusters adjacent to i (candidate moves)

            BestZ = Z[i]                # The default is do nothing
            BestImprove = 0
            V_best = 0
            μ_best = 0
            M_best = 0
            pen_best = 0

            # Set of edges that i is in
            Cv = node2edges[i]

            # Set of hyperedges that node i is in, but don't include i itself
            Cv_list = edgelists[i]

            # Volume of the set currently
            vS = sum(d[Ci])
            dv = d[i]
            @assert(H.D[i] == dv)

            # Check if it's better to move to a nearby cluster, Cj
            for j = 1:length(NC)
                Cj_ind = NC[j]

                # Check improvement for move from cluster i to to cluster j
                if Cj_ind == Ci_ind
                    change = 0
                else
                    ΔV, Δμ, ΔM = increments(V, μ, M, [i], Cj_ind, H.D, Z)
                    voldiff = sum(Ω.ω(Ω.aggregator(p),α)*ΔM[p]*C[p] for p in keys(ΔM))

                    cdiff, Δpen = CutDiff_OneNode(Hyp,w,node2edges,Z,i,Cj_ind,edge2penalty,Ω;α=α) # NOTE: need to check on this
                    change =  cdiff - voldiff


                    ## Alternate way

                    # Change in volume
                    Cj = Clusters[Cj_ind]       # The cluster
                    vJ = sum(d[Cj])
                    Δvol = 0
                    for k = 1:kmax
                        Δvol += bet[k]*((vS-dv)^k + (vJ+dv)^k - vS^k - vJ^k)  # Better if this is smaller
                    end

                    # Change in cut
                    Δcut = 0
                    for eid = 1:length(Cv)
                        # change in cut for edge e when v moves from
                        # cluster Ci to cluster Cj
                        e = Cv[eid]
                        edge_noi = Cv_list[eid]
                        k = elen[e]      # size of the edge
                        we = alp[k]*w[e]
                        if k > 1
                            mc = move_cut(i,Z,edge_noi,Ci_ind,Cj_ind,we)
                        else
                            mc = 0
                        end
                        Δcut += mc
                    end

                    # if abs(cdiff+Δcut) > 1e-8
                    #     println("cut: $i,  $cdiff \t $Δcut")
                    # end
                    #
                    if abs(voldiff-Δvol) > 1e-10 && checkvols
                        println("vol: $i,  old = $voldiff \t new = $Δvol))")
                        # @show volvec
                        # println("vol: $i,  $(abs(voldiff-Δvol))")
                    end

                end

                # Check if this is currently the best possible greedy move
                if change - BestImprove > 1e-8
                    V_best = ΔV
                    μ_best = Δμ
                    M_best = ΔM
                    pen_best = Δpen
                    BestImprove = change
                    BestZ = Cj_ind
                    improving = true
                end
            end

            # Move i to the best new cluster, if it strictly improves modularity
            if BestImprove > 1e-8

                changemade = true
                improving = true    # we have a reason to keep iterating!

                # update increments for volume computation
                V, μ, M = addIncrements(V, μ, M, V_best, μ_best, M_best)

                # update clustering
                ci_old = Z[i]
                Z[i] = BestZ

                # Remove i from its old cluster and add it to its new cluster
                Clusters[ci_old] = setdiff(Clusters[ci_old],i)
                push!(Clusters[BestZ],i)

                # update map from edge ID to penalty at edge
                for eid in keys(pen_best)
                    newpenalty = pen_best[eid]
                    edge2penalty[eid] = newpenalty
                end
            end
        end
    end
    mainloop = time()-tic

    if ~changemade
        if verbose println("No nodes moved clusters") end
    end
    Z, Clusters = renumber(Z,Clusters)
    @show mainloop
    return Z

end

function compute_moddiff(edge2part,Cts,Hyp,w,node2edges,V::Array, μ::Array, M::Dict,I_::T, t::Int64, D::Vector{Int64}, Z::Vector{Int64},C::Dict,Ω::IntensityFunction;α) where T <: Union{Int64, Vector{Int64}}
    """
    NOTE: I_ can now be a Vector{Int64} of nodes, all of which are assumed to belong in the same cluster.
    """
    # increments due to proposal
    ΔV, Δμ, ΔM = increments(V, μ, M, I_, t, D, Z)

    voldiff = sum(Ω.ω(Ω.aggregator(p),α)*ΔM[p]*C[p] for p in keys(M))

    ΔC, Δe2p = CutDiff_Many(Cts, I_, t, Z, Hyp, w, node2edges,edge2part)
    cdiff = sum(ΔC[p]*log(Ω.ω(Ω.aggregator(p), α)) for p in keys(ΔC))
    mdiff = cdiff-voldiff

    return mdiff, ΔV, Δμ, ΔM, ΔC, Δe2p
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
