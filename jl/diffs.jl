
# The purpose of the functions in this module is to streamline the computation of the change in modularity associated with moving a group of nodes, assumed to belong to the same cluster (i.e. a "super node") to a different cluster.

include("cut.jl")
include("vol.jl")
include("dict_ops.jl")

function cutDiff(C, I, t, Z, Hyp, w, node2edges)
    """
    C: Dict{Vector{Int64}, Int64}, the dict of current cut values as would be produced by evalCuts().
        C is keyed by partition vectors p.
        C[p] is the # of edges with label partition p.
        Note: Ω is NOT used in the calculation of C
    I: Vector{Int64}, the list of nodes to move, assumed to be in the same group
    t: Int64, the new group to which to move the nodes I
    Z: Vector{Int64}, the current group assignments
    Hyp: list of nodes in each hyperedge
    w:  weight of the hyperedge (number of times this set of nodes is assigned a hyperedge)
    node2edges: node ID to list of edges the node is in.
    return: ΔC, a Dict{Vector{Int64}, Int64}, where ΔC[p] is the change in C[p] caused by moving the nodes in I to group t
    """

    E_id = Vector{Int64}()  # edge IDs that invole some node from I
    for i in I
        append!(E_id,node2edges[i])
    end
    E = unique(E_id)

    ΔC = Dict(p => 0 for i=1:kmax for j = 1:i for p in keys(C))

    Z_prop = copy(Z)
    Z_prop[I] .= t
    for eid in E_id

        # get nodes in edge, and its weight
        e = Hyp[eid]
        we = w[eid]

        # Updated
        ΔC[partitionize(Z[e])]      -= we

        if haskey(ΔC,partitionize(Z_prop[e]))
            ΔC[partitionize(Z_prop[e])] += we
        else
            ΔC[partitionize(Z_prop[e])] = we
        end

    end
    return(ΔC)
end

function cutDiff_template(C, I, t, Z, H)
    """
    NV -- This is the exact format this function was in before I messed with it.
    C: Dict{Vector{Int64}, Int64}, the dict of current cut values as would be produced by evalCuts().
        C is keyed by partition vectors p.
        C[p] is the # of edges with label partition p.
        Note: Ω is NOT used in the calculation of C
    I: Vector{Int64}, the list of nodes to move, assumed to be in the same group
    t: Int64, the new group to which to move the nodes I
    Z: Vector{Int64}, the current group assignments
    H: the hypergraph, in whatever data format is convenient
    return: ΔC, a Dict{Vector{Int64}, Int64}, where ΔC[p] is the change in C[p] caused by moving the nodes in I to group t
    """

    println("Not implemented. Nate, can you help on this? I think the first line is the main issue")

    # help needed here
    # form list of edges E, where E consists of edges incident to one or more nodes I

    # I think the rest should be ok

    # ΔC = Dict(p => 0 for i=1:kmax for j = 1:i for p in keys(C))

    # Z_prop = copy(Z)
    # Z_prop[I] .= t
    # for e in E
    #     ΔC[partitionize(Z[e])]      -= 1
    #     ΔC[partitionize(Z_prop[e])] += 1
    # end
    # return(ΔC)
end

function volDiff(v, μ, M, I, t, D, Z)
    """
    at the moment, this is just a wrapper for increments() to align with our emerging naming conventions.
    """
    # increments due to proposal
    ΔV, Δμ, ΔM = increments(v, μ, M, I, t, D, Z)
    return ΔV, Δμ, ΔM
end

function modDiff(C, S, I, t, Z, H, v, μ, coefs, Ω, α)
    """
    C: dict of cuts
    S: dict of vols
    I: group of nodes, assumed to be in the same group
    j: proposed label for the nodes in I
    Z: current labeling
    H: hypergraph (probably should be as a sparse matrix)
    v: vector of vols (used in volDiff())
    μ: vector of volume moments (used in volDiff())

    Intended usage (replaces line 487 of hypergraph_louvain.jl)
        ΔQ = modDiff(C, S, I, j, Z, H, v, μ, coefs, Ω, α)
    """

    ΔC          = cutDiff(C, I, t, Z, H)
    ΔV, Δμ, ΔM  = volDiff(S, I, t, Z, H, v, μ)

    # could streamline this with some syntactic sugar if we wanted
    diff = 0.0
    for p in keys(V_prop)
        Op    = Ω(p;α=α,mode="partition")
        diff += log(Op)*ΔC[p] - coefs[p]*ΔM*Op
    end
    return diff
end
