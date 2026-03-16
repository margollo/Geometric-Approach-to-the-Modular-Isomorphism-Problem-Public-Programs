## functions to produce the results for the groups of order 64. The output is a table written in the file 64.tsv which has a similar format as Table 1 in the paper, but wit some additional technical information

# $ julia main.jl

using Oscar
include("../src/minors_strategy.jl") # main functions to apply the minors strategy
include("../test/runtests.jl")       # test functions that the main functions work

"""
An auxiliary function for generating a table.
"""
function record_text(
    i::Int,
    j::Int,
    )::Union{String, Nothing}
    G = small_group(64, i)
    H = small_group(64, j)
    l2 = splitting_layer(G, H, deg=1)
    if (i, j) != (239, 238) # exceptional case
        l4 = splitting_layer(G, H, deg=2)
        l = l4
    else
        l4 = "NA"
        l = splitting_layer(G, H, deg=1)
    end
    Q = practical_geometric_test(G, H, l=l)
    if isnothing(Q)
        return nothing
    else
        return "$i\t$j\t$(small_group_identification(Q))\t$l\t$(l4)\t$(l2)\t$(dimension_modulo(G, l))\t$(dimension_modulo(Q, l))\t$(rank(G))\t$(rank(Q))"
    end
end

const all_pairs = [
    (13, 14),
    (65, 70),
    (105, 104),
    (142, 155),
    (142, 157),
    (155, 157),
    (156, 158),
    (160, 156),
    (160, 158),
    (167, 173),
    (167, 176),
    (176, 173), # about 2.5 hours
    (168, 179),
    (168, 180),
    (180, 179), # about 0.3 hours
    (172, 182),
    (175, 181),
    (239, 238), # the F_4 layer is not available in a reasonable time
]
@assert length(Set(all_pairs)) == 18

# write the table
open("64.tsv", "w") do io
    println(io, "G\tH\tQ\tl\tl(F_4)\tl(F_2)\tdim(Lambda)\tdim(Gamma)\td(G)\td(Q)\tsecond\tGB")
    for (i, j) in all_pairs
        t0 = time_ns()/1e9 # in seconds
        m = round((@allocated r = record_text(i, j))/1e9; digits=2) # in gigabytes
        t1 = time_ns()/1e9
        t = round(t1 - t0; digits=1)
        println(io, "$r\t$t\t$m") 
        flush(io)
    end
end
