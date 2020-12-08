module HypergraphModularity

import Combinatorics
import Distributions
import LinearAlgebra
import Parameters
import Random
import SparseArrays
import StatsBase
import Base
import NLopt
import Optim
import SpecialFunctions
import DelimitedFiles  # 
# eventually move this and the contents of the read_data 
# module over to a separate script, probably in test, 
# as it is not part of the core hypergraph modularity 
# suite

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
include("dyadic.jl")
include("hyper_format.jl")
include("hyperlouvain_helpers.jl")
include("fast_hypergraph_louvain.jl")
include("hypergraph_louvain.jl") # PC: this is deprecated now, right?
include("read_data.jl")          # we should move this out of the package
include("analysis_helpers.jl")   # this too

include("test_funs.jl")

include("AON_hyperlouvain.jl")


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
export read_hypergraph_label_names
export read_hypergraph_labels
export read_hypergraph_edges

export L

export downSampleEdges!
export subhypergraph

export computeDyadicResolutionParameter
export dyadicModularity
export dyadicLogLikelihood

export partitionIntensityFunction
export IntensityFunction
export allOrNothingIntensityFunction
export sumOfExteriorDegreesIntensityFunction

export learnParameters

export AON_Inputs
export SuperNode_PPLouvain

export mutualInformation

export hyperedges
export read_stats_data
export read_cora
export readTemporalData
export subHypergraph
end # module
