## the functions here introduce TableAlgebras, in the sense of the GAP-package ModIsom and defines their arithmetic in Oscar

using Oscar

GAP.Packages.load("ModIsom")

const TableAlgebra = GapObj ### alias for readability
struct TableElem
    tbl::TableAlgebra
    elem::GapObj
end

"""
Return the table element.
"""
function table_elem(tbl::TableAlgebra, A::GapObj)::TableElem
    @assert tbl.dim == length(A)
    @assert all( a in tbl.fld for a in A )
    return TableElem(tbl, A)
end

"""
Return the parent table of `u`.
"""
parent(u::TableElem)::TableAlgebra = u.tbl

"""
Return the underlying element of `u`.
"""
gap(u::TableElem)::GapObj = u.elem

function Base.:(==)(u::TableElem, v::TableElem)::Bool
    return parent(u) == parent(v) && gap(u) == gap(v)
end

function Base.:+(u::TableElem, v::TableElem)::TableElem
    @assert parent(u) == parent(v)
    return table_elem(parent(u), gap(u) + gap(v))
end

Base.:-(u::TableElem)::TableElem = table_elem(parent(u), -gap(u))

Base.:-(u::TableElem, v::TableElem)::TableElem = u + (-v)

"""
This is the usual product. However, the usual symbol is reserved to use `map_word`.
"""
function Base.:%(u::TableElem, v::TableElem)::TableElem
    @assert parent(u) == parent(v)
    return table_elem(parent(u), GAP.Globals.MultByTable(parent(u), gap(u), gap(v)))
end

function Base.show(io::IO, v::TableElem)
    tbl = parent(v)
    for i in 1:tbl.dim
        print(io, string(tbl.pre.exps[i]), "\t")
        print(io, string(tbl.wgs[i]), "\t")
        print(io, string(gap(v)[i]), "\n")
    end
end

"""
This is NOT the usual product.
See u % v for the usual product.
"""
function Base.:*(u::TableElem, v::TableElem)::TableElem
    @assert parent(u) == parent(v)
    return u + v + u % v
end

"""
This is NOT the usual product.
"""
function Base.:^(v::TableElem, n::Int)::TableElem
    @assert n > 0
    return prod(fill(v, n))
end

"""
This is NOT one, but rather zero, which corresponds to one by g -> g - 1.
"""
Base.one(v::TableElem)::TableElem = table_elem(parent(v), 0*gap(v))

"""
Return - v + v^2 - v^3 + ..., which corresponds to inverse by g -> g - 1.
"""
function Base.inv(v::TableElem)::TableElem
    S = [-v]
    lst = -v
    while lst != one(v) # i.e., non-zero element
        lst = (-v) % lst
        push!(S, lst)
    end
    u = sum(S)
    @assert parent(u) == parent(v)
    @assert u * v == one(v)
    @assert v * u == one(v)
    return u
end

function weights(v::TableElem)::Vector{Int}
    tbl = parent(v)
    return Vector{Int}(tbl.wgs)
end

function weights_and_relations(
    v::TableElem,
    gap_to_oscar::MapFromFunc{GapObj, FqMPolyRing},
    )::Vector{Tuple{Int, FqMPolyRingElem}}
    w = weights(v)
    r = gap_to_oscar.(gap(v))
    return [ (wi, ri) for (wi, ri) in zip(w, r) if !is_zero(ri) ]
end
