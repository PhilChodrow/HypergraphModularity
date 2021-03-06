using Revise

using Test

using Pkg; Pkg.activate(".")
using HypergraphModularity
using StatsBase
using Combinatorics

# just testing that these run correctly at the moment
@testset "Ω" begin

    kmax = 10

    # general partitionize-based intensity function
    Ω = partitionIntensityFunction(p->sum(p), kmax)
    @test typeof(Ω) == IntensityFunction

    # all-or-nothing cut
    Ω = allOrNothingIntensityFunction(x->(x[1]+1)^(-x[2]), kmax)
    @test typeof(Ω) == IntensityFunction

    # number of unique elements in z: generalizes all-or-nothing
    Ω = sumOfExteriorDegreesIntensityFunction((x, α)->(x[1]+1)^(-α[x[2]]*x[2]), kmax)
    @test typeof(Ω) == IntensityFunction
end


n = 20
Z = rand(1:5, n)
ϑ = dropdims(ones(1,n) + rand(1,n), dims = 1)
μ = mean(ϑ)
kmax = 4

ω(x, α) = (x[1]+1)^(-α[x[2]]*x[2])

α = repeat([2.0], kmax)

Ω = sumOfExteriorDegreesIntensityFunction(ω, kmax)
H = sampleSBM(Z, ϑ, Ω;α=α, kmax=kmax, kmin = 1)

@testset "HSBM Sampling" begin
    @test length(keys(H.E)) == kmax
end

# test for dyadic log likelihood
# currently broken, need to fix. 

##
@testset "dyadic likelihood" begin
    ωᵢ = 1.0
    ωₒ = 0.1
    L₁ = dyadicLogLikelihood(H, Z, ωᵢ, ωₒ; weighted=false, binary=false, constants = false) 
    L₂ = HypergraphModularity.dyadicLogLikelihoodDirect(H, Z, ωᵢ, ωₒ) 
    @test L₁ ≈ L₂
end
##

##
##







@testset "cut" begin
    cut1 = first_term_eval(H, Z, Ω;α=α)
    Hyp, w = HypergraphModularity.hyperedge_formatting(H)
    cut2 = HypergraphModularity.first_term_v2(Hyp,w,Z,Ω;α=α)
    cut3 = HypergraphModularity.first_term_v3(H, Z, Ω;α=α)
    @test cut1 ≈ cut2
    @test cut2 ≈ cut3
    @test cut1 ≈ cut3

# CutDiff appears to be off, unclear whether related to refactor
    @testset "cutdiff" begin
        Hyp, w = hyperedge_formatting(H)
        node2edges = EdgeMap(H)
        I = rand(1:n)
        J = Z[I] + 1
        a = NaiveCutDiff(H,Z,I,J,Ω;α=α)
        b = HypergraphModularity.NaiveCutDiff2(Hyp, w, node2edges,Z,I,J,Ω;α=α)

        c = CutDiff(Hyp,w,node2edges,Z,I,J,Ω; α=α)

        @test a ≈ b
        # @test a ≈ c
        # @test b ≈ c

        # now we'll test multi-node cutdiff
        # this tests the implementation in diffs.jl.

        to_move = findall(==(3), Z)
        t = 1

        C = HypergraphModularity.evalCuts(H, Z,  Ω)

        ΔC = HypergraphModularity.cutDiff(C, to_move, t, Z, Hyp, w, node2edges, Ω)

        Z_ = copy(Z)
        Z_[to_move] .= t

        C_new = HypergraphModularity.evalCuts(H, Z_,  Ω)

        @test C_new == C + ΔC
    end
end

D = H.D

@testset "vol" begin

    @testset "sums" begin

        all_partitions = collect(Iterators.flatten([(partitions(i,j)) for i =1:kmax for j=1:i]))

        p = [3, 1, 1]

        s1 = HypergraphModularity.evalSumNaive(p, Z, D)
        s2 = HypergraphModularity.evalSumNaive2(p, Z, D)

        @test s1 == s2

        sumsNaive = HypergraphModularity.evalSumsNaive(Z,D,kmax)
        sumsPV    = HypergraphModularity.evalSumsPV(Z,D,kmax)
        sums      = HypergraphModularity.evalSums(Z,D,kmax)[3]

        @test all([sumsNaive[p] == sumsPV[p] for p in all_partitions])
        @test all([sumsNaive[p] == sums[p] for p in all_partitions])
    end

    @testset "voldiff" begin

        ℓ = maximum(Z) # number of total groups

        # let's move all the nodes in group 4 to group 5
        to_move = findall(==(4), Z)
        t = 5
        Z_ = copy(Z)

        # first we'll step through, moving each node one-by-one
        V, μ, M = evalSums(Z, D, kmax; constants=false, bigInt=false);
        for i in to_move
            ΔV, Δμ, ΔM = increments(V, μ, M, i, t, D, Z_)
            V, μ, M = HypergraphModularity.addIncrements(V, μ, M, ΔV, Δμ, ΔM)
            # carry out the change in membership
            Z_[i] = t
        end
        step_wise = M

        # next we'll compute the complete set of increments at once
        V, μ, M = evalSums(Z, D, kmax; constants=false, bigInt=false);
        ΔV, Δμ, ΔM = increments(V, μ, M, to_move, t, D, Z)
        V, μ, M = HypergraphModularity.addIncrements(V, μ, M, ΔV, Δμ, ΔM)
        batch = M
        #
        @test all([step_wise[p] == batch[p] for p in keys(M)])
    end

    @testset "second term" begin

        # aggregation
        ℓ = 0
        bigInt = true
        V, μ, M = evalSums(Z, H, ℓ, bigInt)
        M̂ = HypergraphModularity.aggregateSums(M, Ω)
        @test sum(values(M)) ≈ sum(values(M̂))

        s1 = sum(Ω.ω(Ω.aggregator(p), α)*M[p] for p in keys(M))
        s2 = sum(Ω.ω(p̂, α)*M̂[p̂] for p̂ in keys(M̂))

        @test s1 ≈ s2

        s1 = HypergraphModularity.second_term_eval(H, Z, Ω; α=α)
        s2 = HypergraphModularity.naiveSecondTerm(H, Z, Ω; α=α)

        @test s1 ≈ s2
    end
end

@testset "modularity" begin

    Q_naive = HypergraphModularity.modularityNaive(H, Z, Ω;α=α)
    Q = modularity(H, Z, Ω;α=α)
    @test Q_naive ≈ Q
end

@testset "log likelihood" begin
    Ω = estimateΩEmpirically(H, Z; aggregator = p -> [length(p) == 1, sum(p)], bigNums = true)
    L = HypergraphModularity.logLikelihoodNaive(H, Z, Ω; α=α)
    L̄ = logLikelihood(H, Z, Ω;α=α)
    @test Float64(L) ≈ Float64(sum(L̄))
end







@testset "hyperlouvain" begin
    """
    more a test that the code successfully runs than a meaningful test of correctness
    """
    Z_ = HyperLouvain(H,kmax,Ω,100,true;α=α, verbose=true)
    Z_ = SuperNodeLouvain(H,kmax,Ω,100,true;α=α, verbose=true)
end
