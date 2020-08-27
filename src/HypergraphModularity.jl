module HypergraphModularity

import Combinatorics
import Distributions
import LinearAlgebra
import Parameters
import Random
import SparseArrays
import StatsBase
import Base
import Optim

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
include("fast_hypergraph_louvain.jl")
# include("hypergraph_louvain.jl") # PC: this is deprecated now, right? 
include("read_data.jl")
include("analysis_helpers.jl")

export partitionize
export sampleSBM

export counting_coefficient

export CliqueExpansion
export CliqueExpansionModularity
export HyperLouvain
export SuperNodeLouvain
export SuperNodeStep


export hypergraph

export buildΩ

export evalSums
export evalCuts
export evalConstants
export increments
export addIncremens

export hyperedge_formatting
export EdgeMap
export NaiveCutDiff
export CutDiff

export first_term_eval
export hyperedge_formatting
export first_term_v2
export first_term_v3

export logLikelihood
export modularity
export formObjective

export coordinateAscent

export estimateΩEmpirically

export read_hypergraph_data

export L

export downSampleEdges!

export computeDyadicResolutionParameter
export dyadicModularity
export dyadicLogLikelihood

end # module
