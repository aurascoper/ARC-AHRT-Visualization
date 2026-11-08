# ==============================================================================
# src/structures/StructureIO.jl
# Minimal fail-closed mmCIF reader.
#
# Scope: the categories needed to rebuild a deposited biological assembly and
# audit its symmetry. Anything malformed raises; nothing is reinterpreted or
# guessed. The existing fixed-column PDB reader in M7VisualizationBridge cannot
# represent a 300-chain assembly, and its fallback reinterprets raw bytes rather
# than failing, which is the behaviour this module deliberately does not copy.
# ==============================================================================

"""
    CifBlock

One parsed mmCIF data block. `singles` holds `_category.item -> value` scalars;
`loops` holds `category -> (names, rows)` where `rows` is a vector of vectors of
raw string tokens. Null tokens (`.` and `?`) are preserved verbatim so callers
decide what a missing value means.
"""
struct CifBlock
    name::String
    singles::Dict{String,String}
    loops::Dict{String,Tuple{Vector{String},Vector{Vector{String}}}}
    source_path::String
    source_sha256::String
end

const CIF_NULLS = ("." , "?")

is_cif_null(tok::AbstractString) = tok in CIF_NULLS

# --- tokenizer ---------------------------------------------------------------

"""
    _tokenize_cif_line(line) -> Vector{String}

Split one mmCIF line into tokens, honouring single and double quoting. A quote
only opens a token at a token boundary and only closes when followed by
whitespace or end of line, which is what the CIF 1.1 grammar requires; this is
why `5'-end` inside a bare token does not start a quoted string.
"""
function _tokenize_cif_line(line::AbstractString)
    toks = String[]
    i = firstindex(line)
    n = lastindex(line)
    while i <= n
        c = line[i]
        if isspace(c)
            i = nextind(line, i)
            continue
        end
        if c == '#'              # comment runs to end of line
            break
        end
        if c == '\'' || c == '"'
            q = c
            i = nextind(line, i)
            start = i
            closed = false
            while i <= n
                if line[i] == q && (i == n || isspace(line[nextind(line, i)]))
                    push!(toks, String(line[start:prevind(line, i)]))
                    i = nextind(line, i)
                    closed = true
                    break
                end
                i = nextind(line, i)
            end
            closed || throw(ArgumentError("mmCIF parse error: unterminated $(q) quoted value in line: $(line)"))
        else
            start = i
            while i <= n && !isspace(line[i])
                i = nextind(line, i)
            end
            push!(toks, String(line[start:prevind(line, i)]))
        end
    end
    return toks
end

"""
    read_mmcif(path) -> CifBlock

Parse an mmCIF file. Raises `ArgumentError` on a truncated file, an unterminated
quote or text field, or a `loop_` whose row tokens are not a whole multiple of
its column count. It never pads, truncates or reinterprets a short row.
"""
function read_mmcif(path::AbstractString)
    isfile(path) || throw(ArgumentError("mmCIF file not found: $(path)"))
    raw = read(path)
    isempty(raw) && throw(ArgumentError("mmCIF file is empty: $(path)"))
    digest = bytes2hex(SHA.sha256(raw))
    lines = split(String(raw), '\n')

    blockname = ""
    singles = Dict{String,String}()
    loops = Dict{String,Tuple{Vector{String},Vector{Vector{String}}}}()

    i = 1
    nlines = length(lines)
    pending_key = ""          # a key whose value is on a following line

    while i <= nlines
        line = lines[i]
        stripped = strip(line)

        # multi-line semicolon text field
        if startswith(line, ";")
            buf = IOBuffer()
            print(buf, SubString(line, 2))
            i += 1
            closed = false
            while i <= nlines
                if startswith(lines[i], ";")
                    closed = true
                    i += 1
                    break
                end
                println(buf, lines[i])
                i += 1
            end
            closed || throw(ArgumentError("mmCIF parse error: unterminated ';' text field in $(path)"))
            if !isempty(pending_key)
                singles[pending_key] = strip(String(take!(buf)))
                pending_key = ""
            end
            continue
        end

        if isempty(stripped) || startswith(stripped, "#")
            i += 1
            continue
        end

        if startswith(lowercase(stripped), "data_")
            blockname = String(SubString(stripped, 6))
            i += 1
            continue
        end

        if lowercase(stripped) == "loop_"
            names = String[]
            i += 1
            while i <= nlines
                s = strip(lines[i])
                if isempty(s) || startswith(s, "#")
                    i += 1
                    continue
                end
                startswith(s, "_") || break
                push!(names, String(s))
                i += 1
            end
            isempty(names) && throw(ArgumentError("mmCIF parse error: loop_ with no column names in $(path)"))

            flat = String[]
            while i <= nlines
                s = lines[i]
                st = strip(s)
                if isempty(st) || startswith(st, "#")
                    i += 1
                    # a blank or comment line ends the loop only once rows have started
                    isempty(flat) ? continue : break
                end
                (startswith(st, "_") || lowercase(st) == "loop_" || startswith(lowercase(st), "data_")) && break
                if startswith(s, ";")
                    buf = IOBuffer()
                    print(buf, SubString(s, 2))
                    i += 1
                    closed = false
                    while i <= nlines
                        if startswith(lines[i], ";")
                            closed = true; i += 1; break
                        end
                        println(buf, lines[i]); i += 1
                    end
                    closed || throw(ArgumentError("mmCIF parse error: unterminated ';' field inside loop in $(path)"))
                    push!(flat, strip(String(take!(buf))))
                    continue
                end
                append!(flat, _tokenize_cif_line(s))
                i += 1
            end

            ncol = length(names)
            if !isempty(flat) && length(flat) % ncol != 0
                throw(ArgumentError(
                    "mmCIF parse error: loop for $(first(names)) has $(length(flat)) tokens, " *
                    "not a multiple of its $(ncol) columns (file truncated or malformed): $(path)"))
            end
            rows = Vector{Vector{String}}()
            for r in 1:(length(flat) ÷ ncol)
                push!(rows, flat[((r - 1) * ncol + 1):(r * ncol)])
            end
            cat = String(split(first(names), '.')[1])
            if haskey(loops, cat)
                prev_names, prev_rows = loops[cat]
                prev_names == names || throw(ArgumentError("mmCIF parse error: repeated loop $(cat) with different columns in $(path)"))
                append!(prev_rows, rows)
            else
                loops[cat] = (names, rows)
            end
            continue
        end

        if startswith(stripped, "_")
            toks = _tokenize_cif_line(stripped)
            key = toks[1]
            if length(toks) == 1
                pending_key = key           # value follows on the next line
            else
                singles[key] = join(toks[2:end], " ")
            end
            i += 1
            continue
        end

        # a bare value line completing a key seen on the previous line
        if !isempty(pending_key)
            toks = _tokenize_cif_line(stripped)
            isempty(toks) || (singles[pending_key] = join(toks, " "))
            pending_key = ""
            i += 1
            continue
        end

        i += 1
    end

    isempty(blockname) && throw(ArgumentError("mmCIF parse error: no data_ block header in $(path)"))
    return CifBlock(blockname, singles, loops, String(path), digest)
end

# --- typed accessors ---------------------------------------------------------

"""
    cif_single(block, key; required=false, default="") -> String
"""
function cif_single(b::CifBlock, key::AbstractString; required::Bool=false, default::AbstractString="")
    v = get(b.singles, String(key), nothing)
    if v === nothing || is_cif_null(v)
        required && throw(ArgumentError("mmCIF: required key $(key) absent from $(b.source_path)"))
        return String(default)
    end
    return v
end

"""
    cif_loop(block, category) -> (colindex::Dict{String,Int}, rows)

Returns the column-name-to-index map (keyed by the item name after the dot) and
the raw rows. Raises if the category is absent, so a caller can never silently
treat a missing category as an empty one.
"""
function cif_loop(b::CifBlock, category::AbstractString)
    cat = String(category)
    if !haskey(b.loops, cat)
        # a single-row category is stored as scalars, not a loop; synthesize one row
        prefix = cat * "."
        keys_here = sort!([k for k in keys(b.singles) if startswith(k, prefix)])
        isempty(keys_here) && throw(ArgumentError("mmCIF: category $(cat) absent from $(b.source_path)"))
        idx = Dict{String,Int}()
        row = String[]
        for (j, k) in enumerate(keys_here)
            idx[String(SubString(k, length(prefix) + 1))] = j
            push!(row, b.singles[k])
        end
        return idx, [row]
    end
    names, rows = b.loops[cat]
    idx = Dict{String,Int}()
    for (j, n) in enumerate(names)
        parts = split(n, '.')
        length(parts) == 2 || throw(ArgumentError("mmCIF: malformed loop column name $(n)"))
        idx[String(parts[2])] = j
    end
    return idx, rows
end

has_cif_category(b::CifBlock, category::AbstractString) =
    haskey(b.loops, String(category)) || any(startswith(k, String(category) * ".") for k in keys(b.singles))

function _req(idx::Dict{String,Int}, col::AbstractString, cat::AbstractString, path::AbstractString)
    haskey(idx, String(col)) || throw(ArgumentError("mmCIF: $(cat) is missing required column $(col) in $(path)"))
    return idx[String(col)]
end

# NOTE: `something(tryparse(x), throw(...))` would evaluate the throw eagerly and
# fail on every value, so the check is written out.
function _cifparse(::Type{Float64}, tok::AbstractString, what::AbstractString)
    v = tryparse(Float64, tok)
    v === nothing && throw(ArgumentError("mmCIF: $(what) is not a number: '$(tok)'"))
    return v
end

# --- atoms -------------------------------------------------------------------

"""
    AtomSite

One `_atom_site` record, reduced to what the geometry engine uses.
`asym_id` is `label_asym_id` (the deposited chain id that assembly generation
refers to), `auth_asym_id` the author chain id.
"""
struct AtomSite
    group::String        # ATOM or HETATM
    atom_id::String      # label_atom_id, e.g. "CA"
    comp_id::String      # residue / ligand code
    asym_id::String      # label_asym_id
    auth_asym_id::String
    seq_id::String       # label_seq_id as deposited ("." for ligands)
    element::String
    x::Float64
    y::Float64
    z::Float64
end

"""
    read_atom_sites(block) -> Vector{AtomSite}

Raises if `_atom_site` is absent or empty, or if any coordinate fails to parse.
A structure with no atoms is an error, never an empty result.
"""
function read_atom_sites(b::CifBlock)
    idx, rows = cif_loop(b, "_atom_site")
    isempty(rows) && throw(ArgumentError("mmCIF: _atom_site is empty in $(b.source_path)"))
    p = b.source_path
    ig = _req(idx, "group_PDB", "_atom_site", p)
    ia = _req(idx, "label_atom_id", "_atom_site", p)
    ic = _req(idx, "label_comp_id", "_atom_site", p)
    il = _req(idx, "label_asym_id", "_atom_site", p)
    is = _req(idx, "label_seq_id", "_atom_site", p)
    ie = _req(idx, "type_symbol", "_atom_site", p)
    ix = _req(idx, "Cartn_x", "_atom_site", p)
    iy = _req(idx, "Cartn_y", "_atom_site", p)
    iz = _req(idx, "Cartn_z", "_atom_site", p)
    iauth = get(idx, "auth_asym_id", il)

    out = Vector{AtomSite}(undef, length(rows))
    @inbounds for (k, r) in enumerate(rows)
        out[k] = AtomSite(r[ig], r[ia], r[ic], r[il], r[iauth], r[is], r[ie],
                          _cifparse(Float64, r[ix], "Cartn_x"),
                          _cifparse(Float64, r[iy], "Cartn_y"),
                          _cifparse(Float64, r[iz], "Cartn_z"))
    end
    return out
end

# --- assembly operators ------------------------------------------------------

"""
    SymOp

One `_pdbx_struct_oper_list` entry: a 3x3 matrix and a translation, with the
deposited id and type carried along so the audit can report them.
"""
struct SymOp
    id::String
    type::String
    R::Matrix{Float64}
    t::Vector{Float64}
end

"""
    read_oper_list(block) -> Dict{String,SymOp}
"""
function read_oper_list(b::CifBlock)
    idx, rows = cif_loop(b, "_pdbx_struct_oper_list")
    isempty(rows) && throw(ArgumentError("mmCIF: _pdbx_struct_oper_list is empty in $(b.source_path)"))
    p = b.source_path
    iid = _req(idx, "id", "_pdbx_struct_oper_list", p)
    ity = get(idx, "type", 0)
    ops = Dict{String,SymOp}()
    for r in rows
        R = Matrix{Float64}(undef, 3, 3)
        for a in 1:3, bq in 1:3
            col = _req(idx, "matrix[$(a)][$(bq)]", "_pdbx_struct_oper_list", p)
            R[a, bq] = _cifparse(Float64, r[col], "matrix[$(a)][$(bq)]")
        end
        t = Vector{Float64}(undef, 3)
        for a in 1:3
            col = _req(idx, "vector[$(a)]", "_pdbx_struct_oper_list", p)
            t[a] = _cifparse(Float64, r[col], "vector[$(a)]")
        end
        id = r[iid]
        haskey(ops, id) && throw(ArgumentError("mmCIF: duplicate operator id '$(id)' in $(p)"))
        ops[id] = SymOp(id, ity == 0 ? "" : r[ity], R, t)
    end
    return ops
end

"""
    expand_oper_expression(expr) -> Vector{String}

Expand a `_pdbx_struct_assembly_gen.oper_expression` such as `(1-60)`, `1,2,3`
or `(1-60)(61)`. Products of two groups are returned as `"a*b"` pair keys, which
`build_assembly` composes in that order. Raises on anything it does not
understand rather than silently dropping operators.
"""
function expand_oper_expression(expr::AbstractString)
    s = strip(String(expr))
    isempty(s) && throw(ArgumentError("empty oper_expression"))

    function expand_group(g::AbstractString)
        out = String[]
        for part in split(g, ',')
            part = strip(part)
            isempty(part) && continue
            if occursin('-', part) && !startswith(part, "-")
                lo_s, hi_s = split(part, '-'; limit=2)
                lo = tryparse(Int, strip(lo_s)); hi = tryparse(Int, strip(hi_s))
                (lo === nothing || hi === nothing) && throw(ArgumentError("unparseable range '$(part)' in oper_expression '$(expr)'"))
                lo <= hi || throw(ArgumentError("inverted range '$(part)' in oper_expression '$(expr)'"))
                append!(out, string.(lo:hi))
            else
                push!(out, String(part))
            end
        end
        isempty(out) && throw(ArgumentError("oper_expression group '$(g)' expanded to nothing"))
        return out
    end

    groups = Vector{Vector{String}}()
    if startswith(s, "(")
        i = 1
        while i <= lastindex(s)
            if s[i] == '('
                j = findnext(==(')'), s, i)
                j === nothing && throw(ArgumentError("unbalanced parenthesis in oper_expression '$(expr)'"))
                push!(groups, expand_group(SubString(s, nextind(s, i), prevind(s, j))))
                i = nextind(s, j)
            elseif isspace(s[i])
                i = nextind(s, i)
            else
                throw(ArgumentError("unexpected character '$(s[i])' in oper_expression '$(expr)'"))
            end
        end
    else
        push!(groups, expand_group(s))
    end

    isempty(groups) && throw(ArgumentError("oper_expression '$(expr)' expanded to nothing"))
    length(groups) == 1 && return groups[1]
    length(groups) == 2 && return [string(a, "*", b) for a in groups[1] for b in groups[2]]
    throw(ArgumentError("oper_expression '$(expr)' has $(length(groups)) groups; only 1 or 2 are supported"))
end

"""
    AssemblySpec

One `_pdbx_struct_assembly_gen` row paired with its declared oligomeric count.
"""
struct AssemblySpec
    assembly_id::String
    oper_ids::Vector{String}
    asym_ids::Vector{String}
    declared_oligomeric_count::Int
end

"""
    read_assembly_spec(block; assembly_id="1") -> AssemblySpec
"""
function read_assembly_spec(b::CifBlock; assembly_id::AbstractString="1")
    gidx, grows = cif_loop(b, "_pdbx_struct_assembly_gen")
    p = b.source_path
    ia = _req(gidx, "assembly_id", "_pdbx_struct_assembly_gen", p)
    io = _req(gidx, "oper_expression", "_pdbx_struct_assembly_gen", p)
    il = _req(gidx, "asym_id_list", "_pdbx_struct_assembly_gen", p)

    sel = [r for r in grows if r[ia] == String(assembly_id)]
    isempty(sel) && throw(ArgumentError("mmCIF: no assembly_gen row for assembly '$(assembly_id)' in $(p)"))
    length(sel) == 1 || throw(ArgumentError(
        "mmCIF: assembly '$(assembly_id)' in $(p) is generated by $(length(sel)) rows; " *
        "multi-row assemblies are not supported and must not be silently collapsed"))
    row = sel[1]

    opers = expand_oper_expression(row[io])
    asyms = String[strip(x) for x in split(row[il], ',') if !isempty(strip(x))]
    isempty(asyms) && throw(ArgumentError("mmCIF: empty asym_id_list for assembly '$(assembly_id)' in $(p)"))

    declared = -1
    if has_cif_category(b, "_pdbx_struct_assembly")
        aidx, arows = cif_loop(b, "_pdbx_struct_assembly")
        ja = _req(aidx, "id", "_pdbx_struct_assembly", p)
        if haskey(aidx, "oligomeric_count")
            jc = aidx["oligomeric_count"]
            for r in arows
                if r[ja] == String(assembly_id) && !is_cif_null(r[jc])
                    declared = something(tryparse(Int, r[jc]), -1)
                end
            end
        end
    end
    return AssemblySpec(String(assembly_id), opers, asyms, declared)
end

"""
    em_reconstruction_info(block) -> NamedTuple

Resolution, symmetry type and point symmetry as deposited, or empty strings when
the category is absent (e.g. for a crystal structure). Reported, never assumed.
"""
function em_reconstruction_info(b::CifBlock)
    if !has_cif_category(b, "_em_3d_reconstruction")
        return (resolution = NaN, symmetry_type = "", point_symmetry = "")
    end
    idx, rows = cif_loop(b, "_em_3d_reconstruction")
    isempty(rows) && return (resolution = NaN, symmetry_type = "", point_symmetry = "")
    r = rows[1]
    res = NaN
    if haskey(idx, "resolution") && !is_cif_null(r[idx["resolution"]])
        res = something(tryparse(Float64, r[idx["resolution"]]), NaN)
    end
    st = haskey(idx, "symmetry_type") && !is_cif_null(r[idx["symmetry_type"]]) ? r[idx["symmetry_type"]] : ""
    ps = ""
    if has_cif_category(b, "_em_single_particle_entity")
        sidx, srows = cif_loop(b, "_em_single_particle_entity")
        if !isempty(srows) && haskey(sidx, "point_symmetry") && !is_cif_null(srows[1][sidx["point_symmetry"]])
            ps = srows[1][sidx["point_symmetry"]]
        end
    end
    return (resolution = res, symmetry_type = String(st), point_symmetry = String(ps))
end

struct_title(b::CifBlock) = cif_single(b, "_struct.title")
