include("cut.jl")
include("vol.jl")

function Q(H::hypergraph, Z::Array{Int64, 1}, ℓ::Int64, Ω; constant_terms::Bool = false, bigInt=true)
    if constant_terms
        println("constant terms not yet implemented, returning without constants")
    end
    cut = first_term_eval(H, Z, ℓ, Ω)
    vol = second_term_eval(H, Z, ℓ, Ω, bigInt)
    return cut - vol
end