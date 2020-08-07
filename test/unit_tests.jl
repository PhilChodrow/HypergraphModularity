using Test
# using Combinatorics
# using Random

# include("vol.jl")
# include("cut.jl")
# include("omega.jl")
# include("HSBM.jl")
# include("objectives.jl")
include("test_funs.jl")

using HypergraphModularity

Random.seed!(4321)

@testset "permutations" begin
    edge = [1, 2, 2, 3, 3, 3, 4, 4, 4, 4]
    perm1 = counting_coefficient(edge)

    l = length(edge)
    p = partitionize(edge)
    lfac = factorial(l)
    perm2 = lfac  # adjusting for all permutations
    pe = partitionize(edge)
    for i = 1:length(pe)
        perm2 /= factorial(pe[i])
    end

    @test perm1 == perm2
end

# let's make some simple, fake data
Z = [1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5] # group partition
D = [3, 4, 2, 5, 6, 4, 3, 2, 5, 2, 2]; # degree sequence

r = 3

@testset "vol" begin

    s1, s2 = test_sums([1, 2, 3, 4], Z, D)
    @test s1 == s2
    s1, s2 = test_sums([2, 2, 3, 4, 4, 5], Z, D)
    @test s1 == s2

    all_partitions = collect(Iterators.flatten([(partitions(i,j)) for i =1:r for j=1:i]))


    p = [3, 1, 1]

    s1 = evalSumNaive(p, Z, D)
    s2 = evalSumNaive2(p, Z, D)

    @test s1 == s2

    # test that the naive and product-of-volumes formulae give the same results

    sumsNaive = evalSumsNaive(Z,D,r)
    sumsPV    = evalSumsPV(Z,D,r)
    sums      = evalSums(Z,D,r)[3]

    @test all([sumsNaive[p] == sumsPV[p] for p in all_partitions])
    @test all([sumsNaive[p] == sums[p] for p in all_partitions])

    @testset "increments" begin
        # test for increments in the volume term (the second term of the modularity)

        r = 5

        function testUpdates(Z, D, r, rounds=100, check=true, bigInt=false)

            ℓ = maximum(Z) # number of total groups

            V, μ, M = evalSums(Z, D, r; constants=false, bigInt=bigInt);
            C = evalConstants(r)

            N = Dict()

            for k = 1:rounds

                # random proposal swap
                i = rand(1:length(D)) # node to move
                t = rand(1:ℓ) # new group
                s = Z[i] # old group

                # increments due to proposal
                ΔV, Δμ, ΔM = increments(V, μ, M, i, t, D, Z);

                # new quantities (assumes we accept every proposal)
                V, μ, M = addIncrements(V, μ, M, ΔV, Δμ, ΔM)

                # carry out the change in membership
                Z[i] = t

                # multiply by combinatorial factors to get the sums we actually want
                N = Dict(p => M[p]*C[p] for p in keys(M))
            end
            if check
                V̄, μ̄, N̄ = evalSums(Z, D, r; constants=true)
                return all([N[p] == N[p] for p in keys(M)])
            end
        end;
        @test testUpdates(Z, D, 5, 100, true, true)
    end

    @testset "Batch Increment" begin
        Z = [1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5] # group partition
        D = [3, 4, 2, 5, 6, 4, 3, 2, 5, 2, 2]; # degree sequence

        r = 3

        ℓ = maximum(Z) # number of total groups

        # let's move all the nodes in group 4 to group 5
        to_move = findall(==(4), Z)
        t = 5
        Z_ = copy(Z)

        # first we'll step through, moving each node one-by-one
        V, μ, M = evalSums(Z, D, r; constants=false, bigInt=false);
        for i in to_move
            ΔV, Δμ, ΔM = increments(V, μ, M, i, t, D, Z_)
            V, μ, M = addIncrements(V, μ, M, ΔV, Δμ, ΔM)
            # carry out the change in membership
            Z_[i] = t
        end
        step_wise = M

        # next we'll compute the complete set of increments at once
        V, μ, M = evalSums(Z, D, r; constants=false, bigInt=false);
        ΔV, Δμ, ΔM = increments(V, μ, M, to_move, t, D, Z)
        V, μ, M = addIncrements(V, μ, M, ΔV, Δμ, ΔM)
        batch = M
        #
        @test all([step_wise[p] == batch[p] for p in keys(M)])
    end
end

# sample from a small HSBM
n = 10
Z = rand(1:5, n)
ϑ = dropdims(ones(1,n) + rand(1,n), dims = 1)
μ = mean(ϑ)

Z = rand(1:2, n)
ϑ = dropdims(ones(1,n) + rand(1,n), dims = 1)

# defining group intensity function Ω
μ = mean(ϑ)

ω(p,α) = (10 .*μ*sum(p))^(-sum(p))*prod(p.^α)^(1/(sum(p)*α))
α0 = 1

kmax = 4

Ω = buildΩ(ω, α0, kmax)

H = sampleSBM(Z, ϑ, Ω;α=α0, kmax=kmax, kmin = 1)

# test for incremental updates in the cut term (first term of the modularity). NATE, plug in here

@testset "cutupdate" begin

    Hyp, w = hyperedge_formatting(H)
    node2edges = EdgeMap(H)
    I = rand(1:n)
    J = Z[I] + 1
    a = NaiveCutDiff(H,Z,I,J,Ω;α=α0)
    b = CutDiff(Hyp,w,node2edges,Z,I,J,Ω;α=α0)

    @test a ≈ b

end


@testset "cut" begin

    kmin = 1

    cut1 = first_term_eval(H,Z,Ω;α=α0)

    Hyp, w = hyperedge_formatting(H)

    cut2 = first_term_v2(Hyp,w,Z,Ω;α=α0)

    cut3 = first_term_v3(Z,H,Ω;α=α0)

    @test cut1 ≈ cut2
    @test cut2 ≈ cut3
end


@testset "modularity" begin

    # compute the true LL
    trueLogLik = logLikelihood(H, Z, Ω;α=α0)

    # compute the three terms including modularity and check for near equality with the true likelihood

    Q, K, R = L(H, Z, Ω; α=α0, bigInt=false)

    @test Q + K + R ≈ trueLogLik
end

@testset "parameter objective" begin

    Q = modularity(H, Z, Ω; α=α0, bigInt=false)
    Q_ = parameterEstimateObjective(H, Z, Ω; bigInt=false)

    @test Q_(α0) ≈ -Q
end
