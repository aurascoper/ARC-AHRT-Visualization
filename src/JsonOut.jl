# ==============================================================================
# src/JsonOut.jl
# Minimal dependency-free JSON writer. The package stays stdlib-only so it runs
# in an offline depot, matching SmallModularBioReactors.jl's dependency style.
# NaN and Inf are written as null: JSON has no representation for them, and
# emitting a bare NaN token would produce a file no strict parser can read.
# ==============================================================================
module JsonOut
export json_write

_esc(s::AbstractString) = sprint() do io
    print(io, '"')
    for c in s
        c == '"'  ? print(io, "\\\"") :
        c == '\\' ? print(io, "\\\\") :
        c == '\n' ? print(io, "\\n")  :
        c == '\r' ? print(io, "\\r")  :
        c == '\t' ? print(io, "\\t")  :
        (c < ' ' ? print(io, "\\u", lpad(string(UInt16(c); base=16), 4, '0')) : print(io, c))
    end
    print(io, '"')
end

_emit(io, ::Nothing, ind, lvl) = print(io, "null")
_emit(io, x::Bool, ind, lvl)   = print(io, x ? "true" : "false")
function _emit(io, x::Real, ind, lvl)
    (isnan(x) || isinf(x)) ? print(io, "null") :
        (x isa Integer ? print(io, x) : print(io, repr(float(x))))
end
_emit(io, x::AbstractString, ind, lvl) = print(io, _esc(x))
_emit(io, x::Symbol, ind, lvl) = print(io, _esc(String(x)))

function _emit(io, v::AbstractVector, ind, lvl)
    isempty(v) && return print(io, "[]")
    print(io, "[\n")
    for (i, e) in enumerate(v)
        print(io, ind^(lvl+1)); _emit(io, e, ind, lvl+1)
        print(io, i < length(v) ? ",\n" : "\n")
    end
    print(io, ind^lvl, "]")
end

function _emit(io, d::AbstractDict, ind, lvl)
    isempty(d) && return print(io, "{}")
    ks = sort!(collect(keys(d)); by=string)
    print(io, "{\n")
    for (i, k) in enumerate(ks)
        print(io, ind^(lvl+1), _esc(string(k)), ": ")
        _emit(io, d[k], ind, lvl+1)
        print(io, i < length(ks) ? ",\n" : "\n")
    end
    print(io, ind^lvl, "}")
end

_emit(io, x::Tuple, ind, lvl) = _emit(io, collect(x), ind, lvl)
_emit(io, x, ind, lvl) = print(io, _esc(string(x)))

"Write `obj` as indented JSON to `path`."
function json_write(path::AbstractString, obj; indent::Int=2)
    open(path, "w") do io
        _emit(io, obj, " "^indent, 0)
        print(io, "\n")
    end
    return path
end
end # module
