function main()
    country_labels = Dict{String,Int64}()
    city_labels = Dict{String,Int64}()
    open("locations-TrivagoClickout.txt") do f
        open("city-labels.txt", "w") do g1
            open("country-labels.txt", "w") do g2
                for line in eachline(f)
                    city = strip(line)
                    if !haskey(city_labels, city)
                        city_labels[city] = length(city_labels) + 1
                    end
                    country = strip(split(line, ",")[2])
                    if !haskey(country_labels, country)
                        country_labels[country] = length(country_labels) + 1
                    end
                    write(g1, "$(city_labels[city])\n")
                    write(g2, "$(country_labels[country])\n")                    
                end
            end
        end
    end

    for (labels, name) in [(city_labels, "city-names"),
                           (country_labels, "country-names")]
        info = [(v, k) for (k, v) in labels]
        sort!(info)
        open("$name.txt", "w") do f
            for (v, k) in info
                write(f, "$v $k\n")
            end
        end
    end
end
