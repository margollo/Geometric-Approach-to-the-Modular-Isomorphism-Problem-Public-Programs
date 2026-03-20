# main functions to apply the minors strategy

using Oscar, Logging

#disable_logging(Info)  # uncomment to throw away all @info's and just obtain basic result
GAP.Packages.load("ModIsom")

include("./verify_input.jl")    # quotient_ids, satisfy_rels
include("./table_elem.jl")      # TableElem

"""
Construct the ambient polynomial ring F_p[alpha_{a b}] in the sense of the paper with an isomorphism. The variables are named a_1,a_2,...,a_n, b_1,b_2,...b_n,... Here the number of different letters is the size of a minimal generating set of G and n is the dimension of the algebra Delta(F_pQ)/Delta^{l+1}(F_pQ) where the generators of Delta(F_pG) are mapped to.
"""
function ambient_polynomial_ring(
    G::PcGroup,
    Q::PcGroup,
    l::Int;
    )::Tuple{FqMPolyRing, MapFromFunc{FqMPolyRing, GapObj}}

    # Construct `strs`
    m = rank(G)
    @assert m <= 26 "Too many generators. Lack of letters." # to have nicely named variables we allow at most 26-generated groups
    dims = GAP.Globals.DimensionsLoewyFactors(GapObj(Q))    # corresponding to the theory of Jennings
    @assert 2 <= l + 1 <= length(dims)
    n = sum(dims[2:l + 1])
    strs = [ "$(x)_$(i)" for x in ('a':'z')[1:m] for i in 1:n ]

    # Construct `S`
    p = prime_of_pgroup(G)
    S, _ = polynomial_ring(GF(p), strs, internal_ordering=:degrevlex)
    # From [Cox-Little-O'Shea 2015, p. 58]:
    # "the grevlex ordering is the most efficient for computations"
    oscar_to_gap = Oscar.iso_oscar_gap(S)
    return S, oscar_to_gap
end

"""
Return values for generators under a homomorphism with an isomorphism. In particular, construct the TableAlgebra corresponding to Delta(FQ)/Delta^{l+1}(FQ) and polynomial ring corresponding to G and Q at layer l.
"""
function working_values(
    G::PcGroup,
    Q::PcGroup,
    l::Int;
    )::Tuple{Vector{TableElem}, MapFromFunc{GapObj, FqMPolyRing}}

    @info "Setting `tbl`"
    S, oscar_to_gap = ambient_polynomial_ring(G, Q, l)
    R = codomain(oscar_to_gap)
    RQ = GAP.Globals.GroupRing(R, GapObj(Q))
    tbl = GAP.Globals.PreSetTable(RQ)
    GAP.Globals.WordFillTable(tbl, GapObj(l))
    GAP.Globals.TableOfRadByCollection(tbl)

    @info "Defining `vals`"
    vars = GAP.Globals.IndeterminatesOfPolynomialRing(R)
    vals = [
        let s = 1 + (k - 1)*tbl.dim, t = k*tbl.dim
            table_elem(tbl, vars[s:t])
        end
        for k in 1:rank(G)
    ]
    gap_to_oscar = inv(oscar_to_gap)
    return vals, gap_to_oscar
end

"""
Return the zero relations with the ambient polynomial ring.
"""
function zero_relations(
    G::PcGroup,
    Q::PcGroup,
    l::Int,
    rels::Vector{FPGroupElem};
    )::Tuple{Vector{FqMPolyRingElem}, FqMPolyRing}
    vals, gap_to_oscar = working_values(G, Q, l)
    alg_rels = TableElem[ (@info "Calculating an algebraic relator corresponding to a relator $r ..."; map_word(r, vals)) for r in rels ]
    @info "Finished calculationg algebraic relators"
    for (rel, alg_rel) in zip(rels, alg_rels)
        @info "Weight, number of terms and degree from a relator $rel:"
        for (wi, ri) in  weights_and_relations(alg_rel, gap_to_oscar)
            @info "$wi, $(length(terms(ri))), $(total_degree(ri)), $ri"
        end
    end
    zero_rels = [ gap_to_oscar(zero_rel) for alg_rel in alg_rels for zero_rel in gap(alg_rel) ]
    S = codomain(gap_to_oscar)
    @assert internal_ordering(S) == :degrevlex
    return zero_rels, S
end

"""
Return a vector of relators of generators of the group `G`.
# Examples
```julia-repl
julia> relators_of_presentation(cyclic_group(5))
1-element Vector{FPGroupElem}:
 F1^5
```
"""
function relators_of_presentation(G::PcGroup)::Vector{FPGroupElem}
    K, _ = simplified_fp_group(fp_group(G))
    @assert ngens(K) == rank(G)
    return relators(K)
end

"""
Return quotient groups of `H`, but not of `G`.
```julia-repl
julia> map(describe, distinguishing_quotients(cyclic_group(8), dihedral_group(8)))
2-element Vector{String}:
 "C2 x C2"
 "D8"
```
"""
function distinguishing_quotients(G::PcGroup, H::PcGroup)::Vector{PcGroup}
    return [
        small_group(n, i)
        for (n, i) in  setdiff(quotient_ids(H), quotient_ids(G))
    ]
end

"""
Return the depth of layer.
# Examples
```julia-repl
julia> G = dihedral_group(8); H = quaternion_group(8);
julia> splitting_layer(G, H) # over prime field
2
julia> splitting_layer(G, H, deg=2) # over quadratic field
3
```
"""
function splitting_layer(G::PcGroup, H::PcGroup; deg::Int=1)::Int
    n, i = small_group_identification(G)
    m, j = small_group_identification(H)
    @info "MIPSplitGroupsByAlgebras (deg=$deg) ..."
    res = GAP.evalstr("MIPSplitGroupsByAlgebras([SmallGroup($n, $i), SmallGroup($m, $j)], $deg).splits")
    @assert length(res) == 1 # length(res) == 0 if G and H are isomorphic
    return res[1][1]
end

"""
Return a dictionary containing data for the minors strategy, this is to say all the basic objects we need to apply the strategy.
"""
function minors_data(
    G::PcGroup,
    H::PcGroup;
    l::Int=splitting_layer(G, H),
    Q::PcGroup=H,
    rels::Vector{FPGroupElem}=relators_of_presentation(G),
    )::Dict{
        String,
        Union{
            FqMPolyRing,
            MPolyIdeal{FqMPolyRingElem},
            MatElem{FqMPolyRingElem},
            PcGroup,
            Int,
        }
    }

    # verify
    @info "Verifying inputs ..."
    @assert small_group_identification(Q) in setdiff(quotient_ids(H), quotient_ids(G)) "Wrong choice of quotient group"
    @assert satisfy_rels(rels, G, rank(G)) "Wrong relators for G"
    @assert !satisfy_rels(rels, Q, rank(G)) "Wrong relators for Q"

    # polynomials and ring (ModIsom)
    zero_rels, S = zero_relations(G, Q, l, rels)

    # matrix
    M = matrix(S, rank(G), div(ngens(S), rank(G)), gens(S))[:, 1:rank(Q)]

    # ideal
    I = ideal(zero_rels)

    return Dict(
        "ring" => S,
        "ideal" => I,
        "matrix" => M,
        "quotient" => Q,
        "layer" => l,
    )
end

"""
Return `true` if all minors of the representation matrix are contained in the radical, and `false` otherwise.
# Examples
```julia-repl
julia> G = dihedral_group(8); H = quaternion_group(8);
julia> all_minors_in_radical(G, H)
false
julia> all_minors_in_radical(G, H, l=3)
true
julia> Q = H; F = free_group(2); x, y = gens(F); rels = [ x^2, y^2 ];
julia> all_minors_in_radical(G, H, l=3, Q=Q, rels=rels)
true
```
"""
function all_minors_in_radical(
    G::PcGroup,
    H::PcGroup;
    l::Int=splitting_layer(G, H), # if no other layer is specified we use the smallest possible candidate as calculated by ModIsom
    Q::PcGroup=H,
    rels::Vector{FPGroupElem}=relators_of_presentation(G),
    )::Bool

    # input
    @info "# Inputs"
    @info "G = small_group$(small_group_identification(G)) = $(describe(G))"
    @info "Relators of G = $(relators(G))"
    @info "H = small_group$(small_group_identification(H)) = $(describe(H))"
    @info "Relators of H = $(relators(H))"
    @info "Q = small_group$(small_group_identification(Q)) = $(describe(Q))"
    @info "Relators of Q = $(relators(Q))"
    @info "l = $l"
    @info "rels = $rels"

    dict = minors_data(G, H, l=l, Q=Q, rels=rels)
    S, I, M = dict["ring"], dict["ideal"], dict["matrix"]
    zero_rels = gens(I)

    # data
    @info "# Matrix"
    @info M

    # output. Only here actual algebraic geometry happens
    return all( (@info "Checking whether the radical contains the minor $f ..."; radical_membership(f, I)) for f in minors(M, number_of_columns(M)) )
end

"""
Return a distinguishing quotient `Q` for `G` and `H` or `nothing`.
"""
function practical_geometric_test(
    G::PcGroup,
    H::PcGroup;
    l::Int=splitting_layer(G, H),
    rels::Vector{FPGroupElem}=relators_of_presentation(G),
    )::Union{PcGroup, Nothing}
    for Q in distinguishing_quotients(G, H)
        if all_minors_in_radical(G, H, Q=Q, l=l, rels=rels)
            @info "Splitting quotient = small_group$(small_group_identification(Q)) = $(describe(Q))"
            return Q
        end
    end
    return nothing
end

"""
Return the dimension of Delta(F_pG)/Delta^{l + 1}(F_pG).
# Examples
```julia-repl
julia> dimension_modulo(dihedral_group(8), 3)
6
```
"""
function dimension_modulo(
    G::PcGroup,
    l::Int,
    )::Int
    dims = GAP.Globals.DimensionsLoewyFactors(GapObj(G))
    @assert l + 1 <= length(dims)
    return sum(dims[2:l + 1])
end
