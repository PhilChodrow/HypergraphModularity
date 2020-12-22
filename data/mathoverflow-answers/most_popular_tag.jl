using Random

function main()
    counts = Dict{Int64,Int64}()
    open("orig-node-labels-mathoverflow-answers.txt") do f
        for line in eachline(f)
            tags = [parse(Int64, t) for t in split(line, ",")]
            for t in tags
                curr = get(counts, t, 0)
                counts[t] = curr + 1
            end
        end
    end

    open("node-labels-mathoverflow-answers.txt", "w") do g
        open("orig-node-labels-mathoverflow-answers.txt") do f
            for line in eachline(f)
                tags = [parse(Int64, t) for t in split(line, ",")]
                cnts = [counts[t] for t in tags]
                most_popular_tags = tags[findall(cnts .== maximum(cnts))]
                shuffle!(most_popular_tags)
                write(g, "$(most_popular_tags[1])\n")
            end
        end
    end
end
