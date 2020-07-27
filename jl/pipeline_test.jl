using Optim

include("omega.jl")
include("HSBM.jl")
include("hypergraph_louvain.jl")
include("inference.jl");
include("warmstart.jl")

################################################################################
# CONSTRUCT SYNTHETIC DATA
################################################################################

n = 200
Z = rand(1:2, n)
ϑ = dropdims(ones(1,n) + rand(1,n), dims = 1)

kmax = 3

# function ω(p, α)
#     k = sum(p)
#     return sum(p)/sum((p .* (1:length(p)).^α[k])) / n^(α[kmax+k]*k)
# end
# α0 = [1, 1, 1, .5, .4, .7]


# planted partition, very simple
function ω(p, α)
    k = sum(p)
    length(p) == 1 ? (1.0*n)^(-α[k]) : (1.0*n)^(-α[kmax + k])
end

α0 = [1.0, 0.7, 2.1, 2.0, 1.0, 2.4]

Ω = buildΩ(ω, α0, kmax)

H = sampleSBM(Z, ϑ, Ω; α=α0, kmax=kmax, kmin = 1)

println([(k, length(H.E[k])) for k in 1:kmax])

for k = 1:3
    p = mean([length(partitionize(Z[e])) == 1 for e in keys(H.E[k])])
    println("$(100*p) % of edges are within a single group.")
end

################################################################################
# WARM START
################################################################################

# get an initial Ẑ via clique-expansion + dyadic Louvain

Z = CliqueExpansionModularity(H)
println("Warm start with $(length(unique(Z))) groups.")

# compare to:

Z = SuperNodeLouvain(H,kmax,Ω;α=α0)

################################################################################
# MAIN LOOP
################################################################################

function estimateParameters(H, Z, Ω, α0)
    """
    Currently assumes that the parameter vector has 2*kmax entries, and that the first and second halves of the parameter vector mean meaningfully different things.
    """

    ℓ = maximum(k for k in keys(H.E)) # size of largest hyperedge

    C       = evalCuts(Z,H)
    V, μ, S = evalSums(Z,H,ℓ,true);

    α = copy(α0)
    kmax = length(α)÷2

    res = 0

    function objective(α)
        obj = 0.0
        for p in keys(S)
            Op   = Ω(p; α=α, mode="partition")
            obj += get(C, p, 0)*log(Op) - S[p]*Op
        end
        return -Float64(obj, RoundDown) # sign is for minimization
    end

    function objective(α, a, k)
        α_ = copy(α)
        α_[k] = a[1]
        return objective(α_)
    end

    for i = 1:50
        println(-Optim.minimum(res))
        # optimization in γ
        for k = (kmax+1):(2*kmax)
            res = optimize(a -> objective(α, a, k), -3.0, 3.0) # very slow and simple -- no gradient information
            α[k] = Optim.minimizer(res)[1]

        end
        # optimizationin β
        for k = 1:kmax
            res = optimize(a -> objective(α, a, k), -100, 100)
            α[k] = Optim.minimizer(res)[1]
        end
    end

    return α, -Optim.minimum(res)
end

α, ll = estimateParameters(H, Z, Ω, α0)

println(α0)

println(α)


#
α, ll = estimateParameters(H, Z, Ω, α0)
Z = SuperNodeLouvain(H,kmax,Ω;α=α)
