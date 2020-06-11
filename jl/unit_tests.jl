using Test
using Combinatorics
include("vol.jl")

# let's make some simple, fake data
Z = [1, 1, 2, 2, 3, 3, 4, 4, 4, 5, 5] # group partition
D = [3, 4, 2, 5, 6, 4, 3, 2, 5, 2, 2]; # degree sequence

r = 3


@testset "test sums" begin

    all_partitions = collect(Iterators.flatten([(partitions(i,j)) for i =1:r for j=1:i]))

    # test that the naive and product-of-volumes formulae give the same results

    sumsNaive = evalSumsNaive(Z,D,r)
    sumsPV    = evalSumsPV(Z,D,r)
    sums      = evalSums(Z,D,r)[3]

    @test all([sumsNaive[p] == sumsPV[p] for p in all_partitions])
    @test all([sumsNaive[p] == sums[p] for p in all_partitions])

end

@testset "increments" begin
    # test for increments in the volume term (the second term of the modularity)
    @testset "volume term" begin
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
                return any([N[p] == N[p] for p in keys(M)])
            end     
        end;

        @test testUpdates(Z, D, 5, 100, true, true)
    end

    # test for incremental updates in the cut term (first term of the modularity). NATE, plug in here
    @testset "cut term" begin
        @test false
    end
end