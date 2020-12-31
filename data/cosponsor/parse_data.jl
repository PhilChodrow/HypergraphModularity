using CSV
using DataFrames

function read_congress(fname)
    party_map = Dict{Int64,Int64}()
    name_map = Dict{Int64,String}()
    house = CSV.File(fname) |> DataFrame
    for grp in groupby(house, :id)
        party = unique(grp[!, :party])
        if length(party) != 1
            println("more than one party affiliation")
            @show grp
            error(0)
        end
        party = party[1]
        id = grp[1, :id]
        if party != 100 && party != 200
            if id == 13100
                # James L. Buckley (conservative, map to republican)
                party = 200
            elseif id == 29147
                # Bernie Sanders (independent, map to democratic)
                party = 100
            elseif id == 10802
                # Harry F. Byrd, Jr. (independent democrat, map to democrat)
                party = 100
            elseif id == 94240
                # Jim Jeffords after he left Republicans. This is Probably the
                # most complicated, but he mostly voted with the democrats.
                party = 100
            else
                println("nontraditional party affiliation")
                @show grp
                error(0)
            end
        else
            party_map[id] = party
            name_map[id] = grp[1, :name]
        end
    end
    return party_map, name_map
end

function get_members()
    spm, snm = read_congress("senate.csv")
    hpm, hnm = read_congress("house.csv")
    for key in keys(spm) ∩ keys(hpm)
        if spm[key] != hpm[key]
            println("different parties...")
            @show key
            error(0)
        end
    end

    # merge party maps
    pm = Dict{Int64,Int64}()
    for (key, val) in spm
        pm[key] = val ÷ 100
    end
    for (key, val) in hpm
        pm[key] = val ÷ 100
    end

    # merge name maps
    nm = Dict{Int64,String}()
    for (key, val) in snm
        nm[key] = val
    end
    for (key, val) in hnm
        nm[key] = val
    end
    
    return pm, nm
end

# Get out hyperedges of sponsor + cosponsors, partitioned by bill type
function bill_hedges(pm)
    sponsors = readlines("sponsors.txt")
    all_cosponsors = readlines("cosponsors.txt")
    bill_types = [String(bill[1:2]) for bill in readlines("bills.txt")]
    hedges = Dict{String, Vector{Vector{Int64}}}()
    for bill_type in unique(bill_types)
        hedges[bill_type] = Vector{Vector{Int64}}()
    end
    
    for (sponsor, cosponsors, bill_type) in zip(sponsors, all_cosponsors, bill_types)
        bill_sponsors = [sponsor]
        for cosponsor in split(cosponsors)
            push!(bill_sponsors, cosponsor)
        end
        bill_sponsors = filter(sp -> sp != "NA", bill_sponsors)
        bill_sponsors = map(sp -> parse(Int64, sp), bill_sponsors)
        bill_sponsors = filter(sp -> haskey(pm, sp), bill_sponsors)        
        if length(bill_sponsors) > 1
            push!(hedges[bill_type], bill_sponsors)
        end
    end

    return hedges
end

# Map hypergraph node IDs to 1...N
function mapped_hedges(hedges)
    id_map = Dict{Int64,Int64}()
    inv_id_map = Dict{Int64,Int64}()
    function get_key(v)
        if !haskey(id_map, v)
            val = length(id_map) + 1
            id_map[v] = val
            inv_id_map[val] = v
        end
        return id_map[v]
    end

    new_hedges = []
    for hedge in hedges
        new_hedge = [get_key(v) for v in hedge]
        push!(new_hedges, new_hedge)
    end

    return new_hedges, inv_id_map
end

function write_dataset(hedges, pm, nm, basename)
    new_hedges, inv_id_map = mapped_hedges(hedges)

    dir = string("$(basename)-congress-bills")
    if !isdir(dir)
        mkdir(dir)
    end

    # hyperedge data
    open("$dir/hyperedges-$dir.txt", "w") do f
        for hedge in new_hedges
            write(f, join(hedge, ','))
            write(f, "\n")
        end
    end
            
    # names of the integer labels
    open("$dir/label-names-$dir.txt", "w") do f
        write(f, "Democrat\n")
        write(f, "Republican\n")        
    end
    
    # labels for each node
    open("$dir/node-labels-$dir.txt", "w") do f
        for i = 1:length(inv_id_map)
            party = pm[inv_id_map[i]]
            write(f, "$(party)\n")
        end
    end

    # names of the nodes
    open("$dir/node-names-$dir.txt", "w") do f
        for i = 1:length(inv_id_map)
            name = nm[inv_id_map[i]]
            write(f, "$(name)\n")
        end
    end
end

function main()
    pm, nm = get_members()
    all_hedges = bill_hedges(pm)
    for (bt, hedges) in all_hedges
        # Ignore house ammendments (small)
        if bt != "HZ"
            write_dataset(hedges, pm, nm, bt)
        end
    end
end
