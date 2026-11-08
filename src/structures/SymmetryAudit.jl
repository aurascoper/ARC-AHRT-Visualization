# ==============================================================================
# src/structures/SymmetryAudit.jl
# Build a deposited biological assembly and VERIFY its symmetry rather than
# trusting the deposited label.
#
# This is the "isolation of load-bearing symmetry" step. Every claim here is
# computed from the operator matrices and the resulting coordinates:
#   - each operator is SO(3) to tolerance
#   - the operator set closes under composition (it is a group)
#   - each operator's rotation order is recovered from its trace
#   - the triangulation number is DERIVED from the polymer subunit count,
#     never asserted
# ==============================================================================

"""
    ChainInstance

One generated copy of one asym unit chain: which operator made it, whether it is
a polymer, its atom coordinates (3 x natoms) and its centroid.
"""
struct ChainInstance
    asym_id::String
    oper_key::String
    is_polymer::Bool
    coords::Matrix{Float64}     # 3 x natoms
    centroid::Vector{Float64}
end

"""
    Assembly

A fully generated biological assembly plus the provenance needed to re-derive it.
"""
struct Assembly
    entry_id::String
    title::String
    source_sha256::String
    chains::Vector{ChainInstance}
    n_polymer_chains::Int
    declared_oligomeric_count::Int
    point_symmetry_declared::String
    resolution_angstrom::Float64
end

const POLYMER_GROUPS = ("ATOM",)

"""
    build_assembly(block; assembly_id="1") -> Assembly

Apply the deposited operators to the asym units named by the assembly. A chain
is counted as polymer when its atoms are `ATOM` records; `HETATM`-only chains
(ligands, ions, waters) are generated but excluded from the subunit count, which
is why 6TAP yields 360 chain instances but 300 polymer chains.
"""
function build_assembly(b::CifBlock; assembly_id::AbstractString="1")
    spec = read_assembly_spec(b; assembly_id=assembly_id)
    ops = read_oper_list(b)
    atoms = read_atom_sites(b)
    em = em_reconstruction_info(b)

    by_asym = Dict{String,Vector{AtomSite}}()
    for a in atoms
        push!(get!(by_asym, a.asym_id, AtomSite[]), a)
    end

    function resolve(key::AbstractString)
        if occursin('*', key)
            p = split(key, '*')
            length(p) == 2 || throw(ArgumentError("unsupported operator product '$(key)'"))
            for q in p
                haskey(ops, String(q)) || throw(ArgumentError("operator '$(q)' referenced by assembly but absent from _pdbx_struct_oper_list"))
            end
            o1 = ops[String(p[1])]; o2 = ops[String(p[2])]
            return (o1.R * o2.R, o1.R * o2.t + o1.t)
        end
        haskey(ops, String(key)) || throw(ArgumentError("operator '$(key)' referenced by assembly but absent from _pdbx_struct_oper_list"))
        o = ops[String(key)]
        return (o.R, o.t)
    end

    chains = ChainInstance[]
    for opkey in spec.oper_ids
        R, t = resolve(opkey)
        for asym in spec.asym_ids
            haskey(by_asym, asym) || continue      # asym listed but with no atoms deposited
            sites = by_asym[asym]
            n = length(sites)
            C = Matrix{Float64}(undef, 3, n)
            @inbounds for (j, s) in enumerate(sites)
                C[1, j] = s.x; C[2, j] = s.y; C[3, j] = s.z
            end
            C = R * C .+ t
            ispoly = any(s.group in POLYMER_GROUPS for s in sites)
            push!(chains, ChainInstance(asym, String(opkey), ispoly, C, vec(sum(C, dims=2)) ./ n))
        end
    end
    isempty(chains) && throw(ArgumentError("assembly '$(assembly_id)' generated no chains"))

    npoly = count(c -> c.is_polymer, chains)
    return Assembly(b.name, struct_title(b), b.source_sha256, chains, npoly,
                    spec.declared_oligomeric_count, em.point_symmetry, em.resolution)
end

# --- symmetry verification ---------------------------------------------------

"""
    so3_residual(R) -> (orthonormality_error, det_error)

How far `R` is from a proper rotation. Reported as a number, not a pass flag, so
the audit can show the worst case.
"""
function so3_residual(R::AbstractMatrix{<:Real})
    orth = maximum(abs.(R' * R - Matrix{Float64}(I, 3, 3)))
    return (orth, abs(det(R) - 1.0))
end

"""
    rotation_order(R; maxorder=12, atol=1e-3) -> Int

Recover n from the rotation angle: trace(R) = 1 + 2cos(theta), theta = 2*pi/n.
Returns 1 for the identity and 0 when no order <= maxorder matches, which is
itself a finding rather than an error.
"""
function rotation_order(R::AbstractMatrix{<:Real}; maxorder::Int=12, atol::Float64=1e-3)
    c = clamp((tr(R) - 1.0) / 2.0, -1.0, 1.0)
    theta = acos(c)
    theta < atol && return 1
    # The order is the smallest n with n*theta a whole multiple of 2pi. Matching
    # only theta == 2pi/n would miss non-primitive rotations: a five-fold axis
    # carries both 2pi/5 and 4pi/5, and the latter is still an order-5 element.
    for n in 2:maxorder
        k = n * theta / (2pi)
        abs(k - round(k)) < atol && return n
    end
    return 0
end

"""
    audit_symmetry(block; assembly_id="1", atol=1e-4) -> Dict

Verify, do not trust. Checks that every operator is SO(3), that the set closes
under composition (each product lands on another member within `atol`), and
reports the order spectrum. For a true icosahedral group of order 60 the orders
are 1 x identity, 15 x 2-fold, 20 x 3-fold, 24 x 5-fold.
"""
function audit_symmetry(b::CifBlock; assembly_id::AbstractString="1", atol::Float64=1e-4)
    spec = read_assembly_spec(b; assembly_id=assembly_id)
    ops = read_oper_list(b)
    keys_used = [k for k in spec.oper_ids if !occursin('*', k)]
    mats = [ops[k].R for k in keys_used]
    n = length(mats)

    worst_orth = 0.0; worst_det = 0.0
    for R in mats
        o, d = so3_residual(R)
        worst_orth = max(worst_orth, o); worst_det = max(worst_det, d)
    end

    # closure: every product must match some member
    closes = true
    worst_closure = 0.0
    unmatched = 0
    if n > 1 && n <= 128
        for A in mats, B in mats
            P = A * B
            best = Inf
            for M in mats
                best = min(best, maximum(abs.(P - M)))
                best <= atol && break
            end
            worst_closure = max(worst_closure, best)
            if best > atol
                closes = false; unmatched += 1
            end
        end
    else
        closes = false
    end

    orders = [rotation_order(R) for R in mats]
    spectrum = Dict{Int,Int}()
    for o in orders
        spectrum[o] = get(spectrum, o, 0) + 1
    end

    ico_reference = Dict(1 => 1, 2 => 15, 3 => 20, 5 => 24)
    is_icosahedral_order60 = (n == 60) && closes && (spectrum == ico_reference)

    return Dict(
        "entry_id" => b.name,
        "n_operators" => n,
        "max_orthonormality_error_dimensionless" => worst_orth,
        "max_determinant_error_dimensionless" => worst_det,
        "all_operators_are_SO3" => (worst_orth <= atol && worst_det <= atol),
        "group_closes_under_composition" => closes,
        "max_closure_residual_dimensionless" => worst_closure,
        "unmatched_products_count" => unmatched,
        "rotation_order_spectrum_count_by_order" => Dict(string(k) => v for (k, v) in spectrum),
        "icosahedral_order60_verified" => is_icosahedral_order60,
        "point_symmetry_declared" => em_reconstruction_info(b).point_symmetry,
        "tolerance_dimensionless" => atol,
        "units" => Dict(
            "max_orthonormality_error_dimensionless" => "dimensionless (matrix element)",
            "max_determinant_error_dimensionless" => "dimensionless",
            "max_closure_residual_dimensionless" => "dimensionless (matrix element)",
            "tolerance_dimensionless" => "dimensionless",
            "n_operators" => "count",
            "unmatched_products_count" => "count",
        ),
    )
end

# --- triangulation number, derived not asserted ------------------------------

"Valid Caspar-Klug triangulation numbers T = h^2 + hk + k^2, ascending."
function ck_triangulation_numbers(maxT::Int=49)
    s = Set{Int}()
    for h in 0:12, k in 0:12
        (h == 0 && k == 0) && continue
        T = h^2 + h * k + k^2
        T <= maxT && push!(s, T)
    end
    return sort!(collect(s))
end

"""
    derive_triangulation(n_polymer_chains) -> Dict

N = 60T for a Caspar-Klug shell, so T = N/60 and T must be a valid
triangulation number. A subunit count that divides by 60 into a non-CK value
means the assembly is NOT a single quasi-equivalent lattice, and the result says
so rather than rounding to the nearest plausible T.
"""
function derive_triangulation(n_polymer_chains::Int; declared_oligomeric_count::Int=-1,
                              chains_per_subunit::Int=1)
    valid = ck_triangulation_numbers()
    # Chain count is not subunit count. A cleaved capsid protein is deposited as
    # several chains per subunit (NwV alpha -> beta + gamma), so prefer the
    # deposited oligomeric count when it is present and report any discrepancy
    # rather than silently picking one.
    inferred = chains_per_subunit > 1 ? n_polymer_chains ÷ chains_per_subunit : n_polymer_chains
    n_sub = declared_oligomeric_count > 0 ? declared_oligomeric_count : inferred
    source = declared_oligomeric_count > 0 ? "declared_oligomeric_count" : "polymer_chain_count"
    mismatch = declared_oligomeric_count > 0 && declared_oligomeric_count != n_polymer_chains
    implied_cps = (declared_oligomeric_count > 0 && declared_oligomeric_count > 0 &&
                   n_polymer_chains % declared_oligomeric_count == 0) ?
                   n_polymer_chains ÷ declared_oligomeric_count : -1

    divisible = (n_sub % 60 == 0)
    q = divisible ? n_sub ÷ 60 : -1
    is_ck = divisible && (q in valid)
    decomps = Vector{Vector{Int}}()
    if divisible && !is_ck && q > 0
        for t in valid
            t >= q && continue
            r = q - t
            (r == 1 || r in valid) && push!(decomps, sort([t, r]; rev=true))
        end
        unique!(decomps)
    end
    return Dict(
        "n_polymer_chains_count" => n_polymer_chains,
        "declared_oligomeric_count" => declared_oligomeric_count,
        "subunit_count_used" => n_sub,
        "subunit_count_source" => source,
        "chain_count_differs_from_declared" => mismatch,
        "implied_chains_per_subunit_count" => implied_cps,
        "divisible_by_60" => divisible,
        "N_over_60_dimensionless" => q,
        "is_single_caspar_klug_lattice" => is_ck,
        "derived_T_dimensionless" => is_ck ? q : nothing,
        "candidate_decompositions_T" => decomps,
        "note" => is_ck ?
            "N = 60T with T a valid Caspar-Klug number (T = h^2 + hk + k^2)." :
            "N/60 is not a valid Caspar-Klug triangulation number, so the assembly is not a single quasi-equivalent lattice; see candidate_decompositions_T.",
        "units" => Dict("all counts" => "count", "T" => "dimensionless"),
    )
end
