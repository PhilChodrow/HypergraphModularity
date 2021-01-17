using Plots
pyplot()
using Random

include("hexbinplots_orig.jl")

function shapecoords(hx::HexBinPlots.HexHistogram)
    h,vh = HexBinPlots.make_shapes(hx)
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
    return (xsh, ysh, vh)
end

using Interpolations
function control_point(xi, xj, yi, yj, dist_from_mid)
    xmid = 0.5 * (xi+xj)
    ymid = 0.5 * (yi+yj)
    # get the angle of y relative to x
    theta = atan((yj-yi) / (xj-xi)) + 0.5pi
    # dist = sqrt((xj-xi)^2 + (yj-yi)^2)
    # dist_from_mid = curvature_scalar * 0.5dist
    # now we have polar coords, we can compute the position, adding to the midpoint
    (xmid + dist_from_mid * cos(theta),
     ymid + dist_from_mid * sin(theta))
end

function _add_edge_points!(xlines,ylines,xsi,ysi,xdi,ydi,curve_amount,interplen)
    xpt, ypt = control_point(xsi,xdi,ysi,ydi,curve_amount)
    xpts = [xsi, xpt, xdi]
    ypts = [ysi, ypt, ydi]
    t = range(0, stop=1, length=3)
    A = hcat(xpts, ypts)
    itp = scale(interpolate(A, BSpline(Cubic(Natural(OnGrid())))), t, 1:2)
    tfine = range(0, 1, length=interplen)
    for t in tfine
        push!(xlines, itp(t, 1))
    end
    for t in tfine
        push!(ylines, itp(t, 2))
    end
end
# Mutating version from HexBinPlots
function _xy2counts!(counts::Dict{(Tuple{Int, Int}), Int},
                     x::AbstractArray,y::AbstractArray,
                     xsize::Real,ysize::Real,x0,y0)
    @inbounds for i in eachindex(x)
        h = convert(HexagonOffsetOddR,
                    cube_round(x[i] - x0,y[i] - y0,xsize, ysize))
        idx = (h.q, h.r)
        counts[idx] = 1 + get(counts,idx,0)
    end
    counts
end

using Statistics
using Hexagons

# Generate counts for a graph edge plot histogram
# hhparam = (xsize,ysize,x0,y0) are the hex-histogrma params
function graphedgehist(src,dst,x,y,ptscale,hhparams,curvature = 0.05)
    xlines = similar(x,0)
    ylines = similar(y,0)
    sizehint!(xlines,max(length(src)+1000,10000000)) # we'll process them in batches this big
    sizehint!(ylines,max(length(src)+1000,10000000))
    # this is the counts of the edge vectors...
    edgecounts = Dict{(Tuple{Int, Int}), Int}()
    # now, run through each edges...
    for ei in eachindex(src)
        si = src[ei]
        di = dst[ei]
        xsi = x[si]
        ysi = y[si]
        xdi = x[di]
        ydi = y[di]
        curve_amount = curvature*sign((di-si))
        edgelen = sqrt((xdi-xsi)^2 + (ydi-ysi)^2)
        interplen = max(ceil(Int,2*edgelen/ptscale),2)
        if edgelen == 0
            continue
        end
        _add_edge_points!(xlines, ylines, xsi, ysi, xdi, ydi, curve_amount, interplen)
        if length(xlines) >= 9500000
            # compress all to edgecounts
            println("Compressing counts... $(ei) / $(length(src))")
            _xy2counts!(edgecounts, xlines, ylines, hhparams...)
            resize!(xlines,0)
            resize!(ylines,0)
        end
    end
    # handle final update
    return _xy2counts!(edgecounts, xlines, ylines, hhparams...)
end

function hexbingraphplot(src,dst,x,y;nbins=500)
    xbins = nbins
    ybins = nbins
    # compute params for hex histogram bins
    xmin,xmax = extrema(x)
    ymin,ymax = extrema(y)
    xspan, yspan = xmax - xmin, ymax - ymin
    xsize, ysize = xspan / xbins, yspan / ybins
    x0, y0 = xmin - xspan / 2,ymin - yspan / 2
    hhparams = (xsize, ysize, x0, y0) # save them as a tuple
    edgescale = 6
    hsize = min((xmax-xmin)/(xbins*edgescale),(ymax-ymin)/(ybins*edgescale))
    # need to compute hhparams
    # get coords for edges
    counts = graphedgehist(src,dst,x,y,hsize,hhparams)
    # do this for the nodes...
    #counts = HexBinPlots.xy2counts_(x,y,hhparams...)
    xh,yh,vh = HexBinPlots.counts2xy_(counts,hhparams...)
    hedges = HexBinPlots.HexHistogram(xh,yh,vh,xsize,ysize,false)
    return hedges
end

function generate_nice_plot(xy, src, dst, imagefilename;
                            labels=[],
                            mymarkeralpha=0.1,
                            mymarkersize=5,
                            invpermval=false,
                            colorscale=[],
                            skipedges=false)
    
    if invpermval
        x = xy[:,1]
        y = xy[:,2]
        x = invperm(sortperm(x))
        y = invperm(sortperm(y))
        xy = Float64.(hcat(x,y))
    end

    if !isempty(labels)
	labels = Int.(labels)
	ulabels = unique(labels)
	d = Dict(sort(ulabels) .=> 1:length(ulabels))
	color_palette = distinguishable_colors(length(ulabels))
        color_palette[end] = RGB(1,1,1);
	colorstouse = map(i->color_palette[d[labels[i]]],1:length(labels))
    end
    
    xd = xy[:,1]
    yd = xy[:,2]

    plot()
    if !skipedges
        hedges = hexbingraphplot(src, dst, xd, yd, nbins=900)
        xsh, ysh, axform = shapecoords(hedges)
        avals = axform/maximum(axform)
        avals = map(x -> x <= 1/4 ? sqrt(x) : ((sqrt(1/4)-1/4) + x), avals )
        avals = 1/1.25.*avals
        plot!(xsh, ysh,
              seriestype=:shape,
              alpha=avals,
              fillcolor=:darkblue,
              legend=false,
              linealpha=0,
              framestyle=:none,
              background=nothing)
    end

    #cscale = colorscale ./ maximum(colorscale)
    #rgbvals = get(ColorSchemes.valentine,cscale)
    #ids = findall(colorscale.!=0)
    ids = collect(length(xd):-1:1)
    scatter!(xd[ids], yd[ids],
             color=colorstouse[ids],
             alpha=mymarkeralpha,
             markersize=mymarkersize,
             legend=false,
             colorbar=false,
             markerstrokewidth=0,
             framestyle=:none,
             background=nothing)

    #plot!(dpi=300,size=(800,800))
    #plot!(dpi=200,size=(600,600))
    savefig(imagefilename)
end

function generate_nice_biplot(xy1, xy2, src, dst, imagefilename;
                              labels=[],
                              mymarkeralpha1=0.1, mymarkeralpha2=0.1,                              
                              mymarkersize1=5, mymarkersize2=5,                              
                              incledges=true,
                              inclsecond=true)
    
    xy = [xy1; xy2]
    xd = xy[:,1]
    yd = xy[:,2]

    plot()
    if incledges
        hedges = hexbingraphplot(src, dst, xd, yd, nbins=900)
        xsh, ysh, axform = shapecoords(hedges)
        avals = axform/maximum(axform)
        avals = map(x -> x <= 1/4 ? sqrt(x) : ((sqrt(1/4)-1/4) + x), avals )
        avals = 1/1.25.*avals
        plot!(xsh, ysh,
              seriestype=:shape,
              alpha=avals,
              fillcolor=:darkblue,
              legend=false,
              linealpha=0,
              framestyle=:none,
              background=nothing)
    end

    if inclsecond
        scatter!(xy2[:, 1], xy2[:, 2],
                 color=RGB(1,1,1),
                 alpha=mymarkeralpha2,
                 markersize=mymarkersize2,
                 legend=false,
                 colorbar=false,
                 framestyle=:none,
                 background=nothing)

    end

    function get_colors()
        if isempty(labels)
            return :red
        end
	labels = Int.(labels)
	ulabels = unique(labels)
	d = Dict(sort(ulabels) .=> 1:length(ulabels))
        # skip first in palette
	color_palette = distinguishable_colors(length(ulabels) + 1)
	return map(i->color_palette[d[labels[i]] + 1], 1:length(labels))
    end
    colorstouse = get_colors()
    scatter!(xy1[:, 1], xy1[:, 2],
             color=colorstouse,
             alpha=mymarkeralpha1,
             markersize=mymarkersize1,
             legend=false,
             colorbar=false,
             framestyle=:none,
             background=nothing)


    plot!(dpi=300,size=(800,800))
    #plot!(dpi=200,size=(400,400))
    savefig(imagefilename)
end

