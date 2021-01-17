using BSON
using LinearAlgebra
using SparseArrays

include("glance.jl")

function biplot(dataset, edges, inclsecond=true)
    data = BSON.load("$dataset-xy.bson")
    B = data["B"]
    xy = data["xy"]'
    labels = data["labels"]
    nnodes = length(labels)
    src, dst = findnz(B)[1:2]
    dst .+= nnodes

    # include hyperedges as labeled
    n = length(labels)
    xy1 = xy[1:nnodes, :]
    xy2 = xy[(nnodes + 1):end, :]

    outname = "viz-biplot-$dataset"
    if edges
        outname = "$(outname)-edges"
    end
    generate_nice_biplot(xy1, xy2,
                         src, dst,
                         "$outname.png",
                         mymarkersize1=3, mymarkersize2=1,
                         mymarkeralpha1=0.5, mymarkeralpha2=0.05,
                         incledges=edges,
                         inclsecond=inclsecond,
                         labels=labels,
                         )
end
