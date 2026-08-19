# Tables.jl integration: a vector of results is a table, one row per result.
#
# Broadcasting a solver hands back a `Vector{CarbonateResult}`, so that vector is what people
# want in a DataFrame. Without this, Tables.jl has no interface to use and falls back to
# reflecting over `propertynames` — which advertises the struct fields as well as the computed
# values, purely so tab-completion shows both. That produced a frame carrying `val` twice over,
# once nested and once flattened, alongside columns of `Ks` and `nothing`.

using Tables

"""
`σ_`-prefixed column names for an uncertainty NamedTuple.

`@generated` because the names depend only on the type. Building them per row would intern ~38
symbols for every sample, which costs more than the rest of the row put together.
"""
@generated function _uncertainty_column_names(::Type{<:NamedTuple{names}}) where {names}
    return :($(map(name -> Symbol("σ_", name), names)))
end

"""
One table row: the computed values, and the uncertainties beside them as `σ_pHtot` and friends.

Two methods rather than a branch, so the shape of a row is settled by the result's type
parameter and never tested at runtime. `Ks` and `unit` stay out: `Ks` is a nested structure
rather than a column, and `unit` is constant across a broadcast, so a column repeating "umol"
per sample would be waste.
"""
_row(result::CarbonateResult{V, Nothing}) where {V} = getfield(result, :val)

function _row(result::CarbonateResult{V, E}) where {V, E <: NamedTuple}
    uncertainties = getfield(result, :err)
    named = NamedTuple{_uncertainty_column_names(E)}(values(uncertainties))
    return merge(getfield(result, :val), named)
end

Tables.istable(::Type{<:AbstractVector{<:CarbonateResult}}) = true
Tables.rowaccess(::Type{<:AbstractVector{<:CarbonateResult}}) = true

"""
Give every row the same names, filling what a row does not have with `nothing`.

Results in one vector need not agree on their columns: some may carry uncertainties and some
not, and results of different scopes carry different values. Tables.jl resolves a ragged
rowtable by taking the schema from the first row and quietly ignoring wider ones, so a plain
result sitting in front of an uncertain one would drop every `σ_` column without a word. Filling
keeps the data and costs the affected columns a `Union{Nothing, Float64}` element type.

Names are unioned in first-seen order, so the common case — a plain row ahead of an uncertain
one — still reads values first and uncertainties after.

A broadcast produces a concretely-typed vector, whose rows are uniform by construction, and
takes the early return.
"""
function _fill_ragged_rows(rows::AbstractVector{R}) where {R <: NamedTuple}
    isconcretetype(R) && return rows

    names = Tuple(unique(Iterators.flatten(keys(row) for row in rows)))
    return [NamedTuple{names}(map(name -> get(row, name, nothing), names)) for row in rows]
end

"""
The rows, materialised.

`map` rather than a lazy generator so that the result is a `Vector` of one concrete NamedTuple
type: Tables.jl reads the schema straight off that and builds concretely-typed columns, where an
unknown schema would leave it discovering names row by row. The allocation is the same one
`DataFrame([r.val for r in results])` already makes, so this costs nothing that the idiom it
replaces did not.
"""
Tables.rows(results::AbstractVector{<:CarbonateResult}) = _fill_ragged_rows(map(_row, results))
