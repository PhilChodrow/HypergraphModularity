## Load packages
using Random
using MAT
using StatsBase
using Revise
using Pkg; Pkg.activate(".")
using HypergraphModularity
using SparseArrays
using Statistics
include("../../src/AON_hyperlouvain.jl")
include("../../src/synthetic_hypergraphs.jl")
include("../../src/louvain_utils.jl")
include("../../src/graph_louvain.jl")
## Run experiment

kmin = 2
kmax = 4
davg = 10

lower = 4
v = range(lower,stop=6,length=10)
Nvals = round.(Int64,[10^i for i in v])
s = 1

# Louvain parameters
maxits = 100
randflag = false
verbose = false
clusterpenalty = 0
r_sizes = ones(kmax-kmin+1) # hyperedges of different sizes, equally likely
N = maximum(Nvals)

aris = zeros(length(Nvals),s)
nmis = zeros(length(Nvals),s)
runs = zeros(length(Nvals),s)
cnum = zeros(length(Nvals),s)

aris_refine = zeros(length(Nvals),s)
nmis_refine = zeros(length(Nvals),s)
runs_refine = zeros(length(Nvals),s)
cnum_refine = zeros(length(Nvals),s)

aris_dyadic = zeros(length(Nvals),s)
nmis_dyadic = zeros(length(Nvals),s)
runs_dyadic = zeros(length(Nvals),s)
cnum_dyadic = zeros(length(Nvals),s)

aris_urefine = zeros(length(Nvals),s)
nmis_urefine = zeros(length(Nvals),s)
runs_urefine = zeros(length(Nvals),s)
cnum_urefine = zeros(length(Nvals),s)

aris_udyadic = zeros(length(Nvals),s)
nmis_udyadic = zeros(length(Nvals),s)
runs_udyadic = zeros(length(Nvals),s)
cnum_udyadic = zeros(length(Nvals),s)

println("Louvain Runtime Performance")
print(rpad("Method", 20))
print(rpad("n", 10))
print(rpad("ARI", 15))
print(rpad("NMI", 15))
print(rpad("#Clusters", 15))
println(rpad("Runtime", 15))
println(rpad("",  80, "-"))
# Kmax = maximum(Kvals)

clusterpenalty = 0
for ni = 1:length(Nvals)
    if ni == 1
        println("")
    end
    global n
    n = round(Int64,Nvals[ni])
    m = davg*n
    edgeweights = ones(m)

    K = round(Int64,n/200)

    cluster_sizes=ones(K)
    cluster_prefs=ones(K)

    pvals = [.6, 1/(n^2), 1/(n^2)]

    # Draw s samples and cluster them with the parameters
    for sample = 1:s
        He2n, e2n, elen, deg, ground_truth = GenerateHypergraphAll(n,m,K,pvals,kmin,kmax,cluster_sizes,r_sizes,cluster_prefs)
        n2e = incidence2elist(SparseArrays.sparse(He2n'))
        H = Elist_to_Hypergraph(e2n)
        #cut_weights, vol_weights, learntime, α = learn_alpha_wrapper(H,ground_truth,kmax,n)
        cut_weights, vol_weights, omega, EdgesAndCuts = learn_omega_aon(e2n,ground_truth,kmax,deg,n)
        γ̂ = computeDyadicResolutionParameter(H, ground_truth)
        gu = computeDyadicResolutionParameter(H, ground_truth;weighted = false)

        # Run with hypergraph version
        Zwarm = collect(1:n)  # no warm start
        tic = time()
        Zs = SuperNode_PPLouvain(n2e,e2n,edgeweights,deg,elen,cut_weights,vol_weights,kmax,randflag,maxits,verbose,Zwarm,clusterpenalty)
        toc = time()-tic
        Z = Zs[:,end]
        clusts = maximum(Z)
        ARI = ari(ground_truth,Z)
        NMI = nmi(ground_truth,Z)

        aris[ni,sample] = ARI
        nmis[ni,sample] = NMI
        runs[ni,sample] = toc
        cnum[ni,sample] = clusts

        # Run dyadic version, weighted
        tic =  time()
        Z_dyadic = CliqueExpansionModularity(H,γ̂;maxits = maxits)
        toc = time()-tic
        clusts = maximum(Z_dyadic)
        ARI = ari(ground_truth,Z_dyadic)
        NMI = nmi(ground_truth,Z_dyadic)
        aris_dyadic[ni,sample] = ARI
        nmis_dyadic[ni,sample] = NMI
        runs_dyadic[ni,sample] = toc
        cnum_dyadic[ni,sample] = clusts

        # Run dyadic version, unweighted
        tic =  time()
        Z_udyadic = CliqueExpansionModularity(H,gu;maxits = maxits,weighted = false)
        toc = time()-tic
        clusts = maximum(Z_udyadic)
        ARI = ari(ground_truth,Z_udyadic)
        NMI = nmi(ground_truth,Z_udyadic)
        aris_udyadic[ni,sample] = ARI
        nmis_udyadic[ni,sample] = NMI
        runs_udyadic[ni,sample] = toc
        cnum_udyadic[ni,sample] = clusts

        # Refine weighted dyadic output
        Zwarm = Z_dyadic
        tic = time()
        Zs = SuperNode_PPLouvain(n2e,e2n,edgeweights,deg,elen,cut_weights,vol_weights,kmax,randflag,maxits,verbose,Zwarm,clusterpenalty)
        toc = time()-tic
        Z = Zs[:,end]
        clusts = maximum(Z)
        ARI = ari(ground_truth,Z)
        NMI = nmi(ground_truth,Z)

        aris_refine[ni,sample] = ARI
        nmis_refine[ni,sample] = NMI
        runs_refine[ni,sample] = toc
        cnum_refine[ni,sample] = clusts

        # Refine unweighted dyadic output
        Zwarm = Z_udyadic
        tic = time()
        Zs = SuperNode_PPLouvain(n2e,e2n,edgeweights,deg,elen,cut_weights,vol_weights,kmax,randflag,maxits,verbose,Zwarm,clusterpenalty)
        toc = time()-tic
        Z = Zs[:,end]
        clusts = maximum(Z)
        ARI = ari(ground_truth,Z)
        NMI = nmi(ground_truth,Z)

        aris_urefine[ni,sample] = ARI
        nmis_urefine[ni,sample] = NMI
        runs_urefine[ni,sample] = toc
        cnum_urefine[ni,sample] = clusts
    end

    runmed = median(runs_udyadic[ni,:])
    arimed = median(aris_udyadic[ni,:])
    nmimed = median(nmis_udyadic[ni,:])
    cnummed = median(cnum_udyadic[ni,:])
    print(rpad("UDyadic", 20))
    print(rpad("$n", 10))
    print(rpad("$(round(arimed,digits = 5))", 15))
    print(rpad("$(round(nmimed,digits = 5))", 15))
    print(rpad("$(round(cnummed))", 15))
    println(rpad("$(round(runmed; digits=3))", 10))

    runmed = median(runs_dyadic[ni,:])
    arimed = median(aris_dyadic[ni,:])
    nmimed = median(nmis_dyadic[ni,:])
    cnummed = median(cnum_dyadic[ni,:])
    print(rpad("WDyadic", 20))
    print(rpad("$n", 10))
    print(rpad("$(round(arimed,digits = 5))", 15))
    print(rpad("$(round(nmimed,digits = 5))", 15))
    print(rpad("$(round(cnummed))", 15))
    println(rpad("$(round(runmed; digits=3))", 10))


    print(rpad("AONLouv", 20))
    runmed = median(runs[ni,:])
    arimed = median(aris[ni,:])
    nmimed = median(nmis[ni,:])
    cnummed = median(cnum[ni,:])
    print(rpad("$n", 10))
    print(rpad("$(round(arimed,digits = 5))", 15))
    print(rpad("$(round(nmimed,digits = 5))", 15))
    print(rpad("$(round(cnummed))", 15))
    print(rpad("$(round(runmed; digits=3))", 10))
    println(rpad("$(K)", 10))
    # println(rpad("",  80, "-"))

    print(rpad("U-Ref", 20))
    runmed = median(runs_urefine[ni,:])
    arimed = median(aris_urefine[ni,:])
    nmimed = median(nmis_urefine[ni,:])
    cnummed = median(cnum_urefine[ni,:])
    print(rpad("$n", 10))
    print(rpad("$(round(arimed,digits = 5))", 15))
    print(rpad("$(round(nmimed,digits = 5))", 15))
    print(rpad("$(round(cnummed))", 15))
    print(rpad("$(round(runmed; digits=3))", 10))
    println(rpad("$(K)", 10))

    print(rpad("W-Ref", 20))
    runmed = median(runs_refine[ni,:])
    arimed = median(aris_refine[ni,:])
    nmimed = median(nmis_refine[ni,:])
    cnummed = median(cnum_refine[ni,:])
    print(rpad("$n", 10))
    print(rpad("$(round(arimed,digits = 5))", 15))
    print(rpad("$(round(nmimed,digits = 5))", 15))
    print(rpad("$(round(cnummed))", 15))
    print(rpad("$(round(runmed; digits=3))", 10))
    println(rpad("$(K)", 10))
    println(rpad("",  90, "-"))

    matwrite("Output/Simpler_Parameters.mat",
    Dict("aris"=>aris,"nmis"=>nmis,"runs"=>runs,"cnum"=>cnum,"pvals"=>pvals,"Nvals"=>Nvals,
    "aris_dyadic"=>aris_dyadic,"nmis_dyadic"=>nmis_dyadic,"runs_dyadic"=>runs_dyadic,"cnum_dyadic"=>cnum_dyadic,
    "aris_refine"=>aris_refine,"nmis_refine"=>nmis_refine,"runs_refine"=>runs_refine,"cnum_refine"=>cnum_refine,
    "aris_udyadic"=>aris_udyadic,"nmis_udyadic"=>nmis_udyadic,"runs_udyadic"=>runs_udyadic,"cnum_udyadic"=>cnum_udyadic,
    "aris_urefine"=>aris_urefine,"nmis_urefine"=>nmis_urefine,"runs_urefine"=>runs_urefine,"cnum_urefine"=>cnum_urefine
    ))

end
