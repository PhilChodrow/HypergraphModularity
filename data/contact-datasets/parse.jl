function read_single(filename::String)
    ret = Int64[]
    open(filename) do f
        for line in eachline(f)
            push!(ret, parse(Int64, line))
        end
    end
    return ret
end

function read_node_map(filename::String)
    nm = Dict{Int64,Int64}()
    open(filename) do f
        for (i, line) in enumerate(eachline(f))
            if i == 1; continue; end
            new, orig = [parse(Int64, v) for v in split(line)]
            nm[orig] = new
        end
    end
    return nm
end

function read_labels(filename::String, nm, n::Int64, incl_gender::Bool)
    labels = zeros(Int64, n)
    label_map = Dict{AbstractString,Int64}()

    open(filename) do f
        for line in eachline(f)
            data = split(line)
            orig_id = parse(Int64, data[1])
            # Some nodes do not participate in any edge
            if !haskey(nm, orig_id); continue; end
            node = nm[orig_id]

            # Label is either class or class and gender (note: some genders listed as "UNKNOWN")
            label = data[2]
            if incl_gender
                label = string(label, " ", data[3])
            end
            if !haskey(label_map, label)
                label_map[label] = length(label_map) + 1
            end
            labels[node] = label_map[label]
        end
    end

    return labels, label_map
end

function main()
    for dataset in ["high-school", "primary-school"]
        for incl_gender_in_label in [true, false]
            # Read data
            nverts = read_single("raw-$dataset/contact-$dataset-nverts.txt")
            simplices = read_single("raw-$dataset/contact-$dataset-simplices.txt")
            nm = read_node_map("raw-$dataset/contact-$dataset-nodemap.txt")
            metadata =
                if     dataset == "high-school";    "metadata_2013"
                elseif dataset == "primary-school"; "metadata_primaryschool"
                else                                ""
                end
            labels, label_map = read_labels("raw-$dataset/$metadata.txt", nm, maximum(simplices),
                                            incl_gender_in_label)

            # Get all unique hyperedges
            hedges = Set{Set{Int64}}()
            let curr_ind = 0
	        for nvert in nverts
    	            simplex = simplices[(curr_ind + 1):(curr_ind + nvert)]
                    push!(hedges, Set{Int64}(simplex))
	            curr_ind += nvert
                end
            end
            
            # setup directory
            dir = string("contact-", dataset, "-classes")
            if incl_gender_in_label
                dir = string(dir, "-gender")
            end
            if !isdir(dir)
                mkdir(dir)
            end
            
            # hyperedge data
            open("$dir/hyperedges-$dir.txt", "w") do f
                for hedge in hedges
                    write(f, join(sort(collect(hedge)), ','))
                    write(f, "\n")
                end
            end
            
            # names of the integer labels
            open("$dir/label-names-$dir.txt", "w") do f
                output = sort([(v, k) for (k, v) in label_map])
                for (v, k) in output
                    write(f, "$k\n")
                end
            end
            
            # labels for each node
            open("$dir/node-labels-$dir.txt", "w") do f
                for l in labels
                    write(f, "$l\n")
                end
            end
        end
    end
end
