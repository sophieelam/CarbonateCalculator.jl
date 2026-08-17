# Input validation for the public entry points.
#
# Julia routes any keyword a function does not declare into its `kwargs...` sink, silently.
# Every entry point in this package has such a sink, so `carbon_system(TA=2300, DIC=2000,
# tempc=2.0)` computes at the default 25 °C and returns a plausible wrong number. That is
# the same trap `_reject_retired_arguments` was added for during the conditions refactor;
# this file generalises the guard to every unrecognised name.

"""
    _declared_keywords(entry_point) -> Vector{Symbol}

The keyword names `entry_point` actually declares, read off its method signature.

Used only to suggest a correction once a call has already failed, so it is never on the hot
path, and reading the signature rather than maintaining a parallel list means the two cannot
drift apart. The `kwargs...` sink itself appears as a name ending in `...` and is dropped.
"""
function _declared_keywords(entry_point)
    names = Symbol[]
    for method in methods(entry_point)
        append!(names, Base.kwarg_decl(method))
    end
    filter!(name -> !endswith(String(name), "..."), names)
    return sort!(unique!(names))
end

"""
Collapse a keyword name to the form that survives the mistakes people actually make —
dropped underscores and inconsistent case. `tempc`, `TempC` and `temp_c` all reduce to
`tempc`, which is enough to catch the common typo without pulling in an edit-distance
implementation for the sake of it.
"""
_normalise_keyword(name::Symbol) = lowercase(replace(String(name), "_" => ""))

"Declared keywords whose normalised form matches `name`, for a 'did you mean' hint."
function _similar_keywords(name::Symbol, entry_point)
    target = _normalise_keyword(name)
    return [k for k in _declared_keywords(entry_point) if _normalise_keyword(k) == target]
end

"""
    _reject_unknown_arguments(kwargs, entry_point)

Reject keywords the calculation does not understand.

`kwargs` is the entry point's keyword sink. Julia routes only *unmatched* keywords there, so
anything present is by definition a name `entry_point` does not accept — the check needs no
list of valid names, and cannot fall out of step with the signature.

Retired condition names are handled first so they keep their more specific message, which
names the replacement.
"""
function _reject_unknown_arguments(kwargs, entry_point)
    _reject_retired_arguments(kwargs)
    isempty(kwargs) && return nothing

    lines = map(collect(keys(kwargs))) do name
        suggestions = _similar_keywords(name, entry_point)
        isempty(suggestions) ? "$name" : "$name (did you mean $(join(suggestions, " or "))?)"
    end

    throw(ArgumentError(
        "unrecognised argument(s) to $(nameof(entry_point)):\n  " * join(lines, "\n  ") *
        "\nUnrecognised keywords used to be absorbed silently, which meant a typo returned " *
        "a plausible result computed at the defaults."
    ))
end
