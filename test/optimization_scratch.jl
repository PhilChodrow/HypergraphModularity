# This simple script illustrates the use of optimization methods to infer values for
# α in a parameterized form of Ω. There are several major issues with the current state of the code.
# - First, we are VERY slow. This could potentially be addressed in part by incorporating gradient information explicitly.
# - Second, the optimization landscape appears to be somewhat rough, and the convergence appears to be very strongly sensitive to initial conditions.
# - Third, the gradients in α[2] (in the current code as written) are likely MUCH steeper than the gradients in α[1], resulting in a very badly conditioned problem that might display poor behavior. Need to find some way to address this, or use a different kind of objective function.

# - On the other hand, the directional results are what we would expect: α[1] is positive, and α[2] > 1.

using StatsBase
using Combinatorics

include("omega.jl")
include("HSBM.jl")
include("read_data.jl")
include("inference.jl")

# cd("hypergraph_modularities_code")

kmax = 10

H, Z = read_hypergraph_data("contact-primary-school",kmax)

n = length(H.D)

# parametric form for Ω

# Let's consider a very simple parameterization of Ω. Much more complex and presumably better ones are possible.

# ω(p,α) = prod(p.^α[1])^(1/(α[1]*length(p))) / (n^(α[2]*sum(p)))

# here's a different one: for each size, we have a different parameter giving dependence on fraction in largest group.

function ω(p, α)
        k = sum(p)
        return (p[1] / k)^α[k]  / n^(α[kmax + k]*k)
end

# α0 = [.6, -1.4] # initial value of α, needed to build Ω and to perform the optimization below.

α0 = repeat([3.0], 2*kmax)

Ω = buildΩ(ω, α0, kmax)

# need more sophisticated control of the optimization, I think.
# basically promising ,but work is needed.
res = estimateΩParametrically(H, Z, Ω, α0,
                              Optim.Options(f_tol      = 1e-10,
                                            iterations = 1000,
                                            show_trace = true))
println(res)

res
