module HexBinPlots
# Copyright (c) 2017: Ralph Smith
# portions derived from Gadfly.jl, StatPlots.jl, StatsBase.jl
# Copyright (c) 2012-2017: the Gadfly.jl developers
# Copyright (c) 2012-2017: the StatPlots.jl developers
# Copyright (c) 2012-2017: the StatsBase.jl developers
# SPDX-License-Identifier: MIT-Expat
# License-Filename: LICENSE.md
export hexbins
using Compat
#using Reexport
#@reexport using Plots
using Hexagons
using StatsBase
import StatsBase.fit
using Plots
# superficially similar to StatsBase.Histogram API
struct HexHistogram{T <: Real,S <: Real}
    xc::Vector{T}
    yc::Vector{T}
    weights::Vector{S}
    xsize::T
    ysize::T
    isdensity::Bool
end
function fit(::Type{HexHistogram},x::AbstractVector,y::AbstractVector,
             bins::Union{NTuple{1,Int},NTuple{2,Int},Int};
             xyweights=nothing,
             density::Bool=false)
    if length(bins) == 2
        xbins, ybins = bins
    else
        xbins, ybins = (bins...,bins...)
    end
    xmin, xmax = extrema(x)
    ymin, ymax = extrema(y)
    xspan, yspan = xmax - xmin, ymax - ymin
    xsize, ysize = xspan / xbins, yspan / ybins
    fit(HexHistogram,x,y,xsize,ysize,boundingbox=[xmin,xmax,ymin,ymax],
        density=density,xyweights=xyweights)
end
function xy2counts_(x::AbstractArray,y::AbstractArray,
                    xsize::Real,ysize::Real,x0,y0)
    counts = Dict{(Tuple{Int, Int}), Int}()
    @inbounds for i in eachindex(x)
        h = convert(HexagonOffsetOddR,
                    cube_round(x[i] - x0,y[i] - y0,xsize, ysize))
        idx = (h.q, h.r)
        counts[idx] = 1 + get(counts,idx,0)
    end
    counts
end
function xy2counts_(x::AbstractArray,y::AbstractArray,wv::AbstractVector{T},
                    xsize::Real,ysize::Real,x0,y0) where{T}
    counts = Dict{(Tuple{Int, Int}), T}()
    zc = zero(T)
    @inbounds for i in eachindex(x)
        h = convert(HexagonOffsetOddR,
                    cube_round(x[i] - x0,y[i] - y0,xsize, ysize))
        idx = (h.q, h.r)
        counts[idx] = wv[i] + get(counts,idx,zc)
    end
    counts
end
function counts2xy_(counts::Dict{S,T},xsize, ysize, x0, y0) where {S,T}
    nhex = length(counts)
    xh = zeros(nhex)
    yh = zeros(nhex)
    vh = zeros(T,nhex)
    k = 0
    for (idx, cnt) in counts
        k += 1
        xx,yy = Hexagons.center(HexagonOffsetOddR(idx[1], idx[2]),
                                xsize, ysize, x0, y0)
        xh[k] = xx
        yh[k] = yy
        vh[k] = cnt
    end
    xh,yh,vh
end
function fit(::Type{HexHistogram},x::AbstractVector,y::AbstractVector,
             xsize, ysize; boundingbox=[], density::Bool=false,
             xyweights::Union{Nothing,AbstractWeights}=nothing)
    (length(x) == length(y)) || throw(
        ArgumentError("data vectors must be commensurate"))
    (xyweights == nothing ) || (length(xyweights) == length(x)) || throw(
        ArgumentError("data and weight vectors must be commensurate"))
    if isempty(boundingbox)
        xmin, xmax = extrema(x)
        ymin, ymax = extrema(y)
    else
        xmin, xmax, ymin, ymax = (boundingbox...,)
    end
    xspan, yspan = xmax - xmin, ymax - ymin
    x0, y0 = xmin - xspan / 2,ymin - yspan / 2
    if xyweights == nothing
        counts = xy2counts_(x,y,xsize,ysize,x0,y0)
    else
        counts = xy2counts_(x,y,xyweights.values,xsize,ysize,x0,y0)
    end
    xh,yh,vh = counts2xy_(counts,xsize, ysize, x0, y0)
    if density
        binarea = sqrt(27)*xsize*ysize/2
        vh = 1/(sum(vh)*binarea) * vh
    end
    HexHistogram(xh,yh,vh,xsize,ysize,density)
end
function make_shapes(h::HexHistogram)
    nhex = length(h.xc)
    hexes = Vector{Shape}(undef, nhex)
    for k in eachindex(h.xc)
        hp = hexpoints(h.xc[k],h.yc[k],h.xsize,h.ysize)
        hexes[k] = Shape([p[1] for p in hp],[p[2] for p in hp])
    end
    hexes,h.weights
end
"""
    hexbins(x,y; density=false, weights=ones(), bins=20,...) --> plot object
Construct and plot a hexbin histogram of two vectors `x`,`y` of equal length.
# Important keyword options
* `bins`: an integer or 2-tuple of integers determining grid resolution.
* `density`: boolean indicating whether result should be normalized as a density.
* `weights`: `AbstractWeights` denoting contribution of each sample pair to its bin.
* Plot attributes, notably `colorbar`, `color_palette`.
"""
function hexbins end
@shorthands hexbins
@recipe function f(::Type{Val{:hexbins}}, plt::AbstractPlot; density = false,
                   weights=nothing)
    x, y = plotattributes[:x], plotattributes[:y]
    bns = get(plotattributes, :bins, 20)
    hx = fit(HexHistogram,x,y,bns,density=density,xyweights=weights)
    legend := false
    h,vh = make_shapes(hx)
    vmax = maximum(vh)
    xsh = Vector{Float64}()
    ysh = Vector{Float64}()
    for k in eachindex(h)
        append!(xsh,h[k].x)
        push!(xsh,h[k].x[1])
        push!(xsh,NaN)
        append!(ysh,h[k].y)
        push!(ysh,h[k].y[1])
        push!(ysh,NaN)
    end
    if haskey(plotattributes,:color_palette)
        fc = Plots.cgrad(plotattributes[:color_palette])
    else
        fc = Plots.cgrad()
    end
    cb = get(plotattributes,:colorbar,false)
    wantcb = (cb != :none) & (cb != false)
    if wantcb && (cb != true) && (cb != :right) && (cb != :best)
        warn("overriding unimplemented colorbar specification")
    end
    if wantcb
        xr = 1:3
        yr = linspace(0,vmax,256)
        ramp = ones(3)' .* yr
        layout --> @layout [ ctrplot{0.95w,1h} cbar ]
    end
    @series begin
        seriestype := :shape
        subplot := 1
        background_color_inside := fc[0.0]
        cols = map(x->fc[clamp(x/vmax,0,1)],vh)
        fillcolor := cols
        x := xsh
        y := ysh
    end
    if wantcb
        @series begin
            seriestype := :heatmap
            # left_margin --> 1mm
            colorbar := false
            subplot := 2
            xticks := false
            ylims := extrema(yr)
            ymirror := true
            if haskey(plotattributes,:color_palette)
                fillcolor := Plots.cgrad(plotattributes[:color_palette])
            end
            x := xr
            y := yr
            z := Surface(ramp)
        end
    end
end
end