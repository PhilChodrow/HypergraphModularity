using Pkg; Pkg.activate(".")
using HypergraphModularity
using StatsBase
using DataFrames
using SparseArrays
using Clustering
using CSV

function recoveryExperiment(; data, n_rounds = 2, core = 0, Γ₀ = 100.0, γ₀ = 1.0, kmax = 20, save_as = "recovery_throughput/experiment.csv")
    
    # collector
    DF = DataFrame()

    H, Z = read_hypergraph_data(data, kmax)
    
    H, Z = kcore(H, Z, core)
    n = length(H.D)
    
    kmin = minimum([k for k in keys(H.E) if length(H.E[k]) > 0])
    
    println("$(core)-core of $(data) has $n nodes")
    
    for k ∈ reverse(kmin:kmax)
        # remove edges of size larger than k
        removeEdges!(H; remove = collect((k+1):(kmax+1)))
        
        m = sum(length(H.E[k]) for k in keys(H.E))
        
        # polyadic warmstart
        D = big(sum(H.D))

        ωᵢ = 2.0
        ωₒ = 0.6
        function ω(p, α)
            num = p[1] == 1 ? ωᵢ : ωₒ
            denom = D^sum(p[2])
            return Γ₀*num / denom
        end
        
        Ω̂ = allOrNothingIntensityFunction(ω, maximum(keys(H.E)))
        Ẑ_w = SuperNode_PPLouvain(H, Ω̂; α = 0, verbose = false, Z0 = collect(1:n));
        
        # polyadic experiments
        for i ∈ 1:n_rounds 
            Ẑ, t = @timed SuperNode_PPLouvain(H, Ω̂; α = 0, verbose = false, scan_order = "random")
            Ω̄ = estimateΩEmpirically(H, Ẑ; aggregator = p -> [length(p) == 1, sum(p)])

            Q = modularity(H, Ẑ, Ω̂; α = 0)
            ℓ = length(unique(Ẑ))    
            ARI = randindex(Z, Ẑ)[1]
            df = DataFrame(data = data, kmax = k, ℓ = ℓ, Q = Q, t = t, ARI = ARI, method = "Polyadic", n = n, core = core, round = i, thread = Threads.threadid(), m = m)
            DF = vcat(DF, df)
        end
        
        
        # construct dyadic graph
        H̄ = projectedGraph(H)  
        
        # dyadic warmstart
        D̄ = big(sum(H̄.D))
        
#         function ω_d(p, α)
#             num = p[1] == 1 ? ωᵢ : ωₒ
#             denom = D̄^sum(p[2])
#             return γ₀*num / denom
#         end
        
#         Ω̄ = allOrNothingIntensityFunction(ω_d, 2)
#         Z̄_w = SuperNode_PPLouvain(H, Ω̄; α = 0, verbose = false, Z0 = collect(1:n));
              
        best_Z̄ = zero(Z)
        best_Q = -Inf
        
        Ω̄ = estimateΩEmpirically(H̄, Ẑ_w; aggregator = p -> [length(p) == 1, sum(p)])
        for i ∈ 1:n_rounds
            Z̄, t = @timed SuperNode_PPLouvain(H̄, Ω̄; α = 0, verbose = false)
            Ω̄ = estimateΩEmpirically(H̄, Z̄; aggregator = p -> [length(p) == 1, sum(p)])
            Q = modularity(H̄, Z̄, Ω̄; α = 0)
            ℓ = length(unique(Z̄))
            ARI = randindex(Z, Z̄)[1]
            df = DataFrame(data = data, kmax = k, ℓ = ℓ, Q = Q, method = "Dyadic", t = t, ARI = ARI, n = n, core = core, round = i, thread = Threads.threadid(), m = m)
            DF = vcat(DF, df)
            
                        
            if Q > best_Q
                best_Q = Q
                best_Z̄ = Z̄
            end   
        end
        
        # dyadic with weighted projection

        γ = computeDyadicResolutionParameter(H, Ẑ_w; mode = "γ", weighted = true)
        
        for i ∈ 1:n_rounds
            Z̄, t = @timed CliqueExpansionModularity(H, γ; weighted = true)
            γ = computeDyadicResolutionParameter(H, Z̄; mode = "γ", weighted = true)
            Ω̄ = estimateΩEmpirically(H̄, Z̄; aggregator = p -> [length(p) == 1, sum(p)])
            Q = modularity(H̄, Z̄, Ω̄; α = 0)
            ℓ = length(unique(Z̄))
            ARI = randindex(Z, Z̄)[1]
            df = DataFrame(data = data, kmax = k, ℓ = ℓ, Q = Q, method = "Dyadic (weighted)", t = t, ARI = ARI, n = n, core = core, round = i, thread = Threads.threadid(), m = m)
            DF = vcat(DF, df)           
        end
        
        # refinement
        Ẑ = best_Z̄    
        for i ∈ 1:n_rounds
            Ω̂ = estimateΩEmpirically(H, Ẑ; aggregator = p -> [length(p) == 1, sum(p)])
            Ẑ, t = @timed SuperNode_PPLouvain(H, Ω̂; Z0 = Ẑ, α = 0, scan_order = "random", verbose = false)
            Q = modularity(H, Ẑ, Ω̂; α = 0)
            ℓ = length(unique(Ẑ))    
            ARI = randindex(Z, Ẑ)[1]
            df = DataFrame(data = data, kmax = k, ℓ = ℓ, Q = Q, t = t, ARI = ARI, method = "Polyadic Refinement", n = n, core = core, round = i, thread = Threads.threadid(), m = m)
            DF = vcat(DF, df)
        end  
        path = "fig/recovery_throughput/$(k)_$(save_as).csv"
        CSV.write(path, DF)
    end
    
    
    return DF
end

## Main experiment
n_rounds = 20
kmax = 25

control = [
    Dict(:data => "TrivagoClickout-fix",   :core => 2,  :n_rounds => n_rounds, :Γ₀ => 1000.0, :γ₀ => 100000.0, :kmax => kmax, :save_as => "TrivagoClickout-fix-2-core"),
    Dict(:data => "TrivagoClickout",   :core => 2,  :n_rounds => n_rounds, :Γ₀ => 1000.0, :γ₀ => 100000.0, :kmax => kmax, :save_as => "TrivagoClickout-2-core")
];

Threads.@threads for d in control
    recoveryExperiment(;d...)
end