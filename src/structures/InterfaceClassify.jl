# ==============================================================================
# src/structures/InterfaceClassify.jl
# Give each contact-graph edge class a STRUCTURAL identity rather than a bare
# centroid separation.
#
# For every edge class we report:
#   - multiplicity, and multiplicity / 60, which is how many such contacts exist
#     per icosahedral asymmetric unit
#   - the chain pair (same deposited chain id, or a cross pair)
#   - the symmetry relation: which rotation order maps one partner onto the
#     other, recovered by finding the operator whose action carries subunit i
#     closest to subunit j
#
# The rotation order is what names the interface. An edge whose partners are
# related by a five-fold operator sits on a pentamer; by a three-fold, on a
# trimer axis; by a two-fold, on a dimer axis. Nothing here is hand assigned.
# ==============================================================================

"""
    classify_interfaces_structurally(block; cutoff_angstrom, assembly_id="1", ndigits=0)

Returns one row per edge class with its multiplicity, chain pairing and the
rotation order of the operator relating the two partners.
"""
function classify_interfaces_structurally(b::CifBlock; cutoff_angstrom::Real,
                                          assembly_id::AbstractString="1", ndigits::Int=0)
    asm = build_assembly(b; assembly_id=assembly_id)
    ops = read_oper_list(b)
    poly = [c for c in asm.chains if c.is_polymer]
    cen = [c.centroid for c in poly]
    n = length(cen)

    # order of each operator, by its deposited id
    oporder = Dict(k => rotation_order(o.R) for (k, o) in ops)

    nv, edges = subunit_contact_graph(cen, cutoff_angstrom)

    # global centroid, so we can test a partner against each operator's action
    C = reduce(hcat, cen)
    mu = vec(sum(C, dims=2)) ./ n

    "rotation order of the operator carrying subunit i closest onto subunit j"
    function relating_order(i::Int, j::Int)
        pi_ = cen[i] .- mu
        pj_ = cen[j] .- mu
        best_k = ""; best_d = Inf
        for (k, o) in ops
            d = norm(o.R * pi_ .- pj_)
            if d < best_d
                best_d = d; best_k = k
            end
        end
        return get(oporder, best_k, 0), best_d, best_k
    end

    groups = Dict{String,Vector{Int}}()
    for (e, (u, v)) in enumerate(edges)
        d = norm(cen[u] .- cen[v])
        push!(get!(groups, string(round(d; digits=ndigits)), Int[]), e)
    end

    rows = Dict{String,Any}[]
    for (key, idxs) in groups
        u, v = edges[idxs[1]]
        ord, resid, opk = relating_order(u, v)
        chains = sort(unique([poly[u].asym_id, poly[v].asym_id]))
        # how many partners of this class each subunit has
        deg = zeros(Int, n)
        for e in idxs
            a, bq = edges[e]; deg[a] += 1; deg[bq] += 1
        end
        sep = parse(Float64, key)
        push!(rows, Dict{String,Any}(
            "separation_angstrom" => sep,
            "separation_nm" => sep / 10,
            "multiplicity_count" => length(idxs),
            "copies_per_asymmetric_unit" => length(idxs) / 60,
            "chain_pair" => join(chains, "-"),
            "same_chain_type" => length(chains) == 1,
            "relating_rotation_order" => ord,
            "relating_operator_id" => opk,
            "superposition_residual_angstrom" => resid,
            "mean_partners_per_subunit" => sum(deg) / n,
            # The operator test is only meaningful when a rotation actually
            # carries one partner onto the other, which requires them to be
            # copies of the SAME deposited chain. For a cross-chain contact no
            # operator can superpose them, the search returns a large residual,
            # and the recovered "order" is an artifact of the fallback. Those
            # classes are named by their chain pairing and left for a human to
            # assign to quasi-equivalent positions.
            "axis_name" => resid > 2.0 ?
                "cross-subunit contact, not related by a rotation (NEEDS NAMING)" :
                (ord == 5 ? "five-fold (pentamer)" :
                 ord == 3 ? "three-fold (trimer)" :
                 ord == 2 ? "two-fold (dimer)" :
                 ord == 1 ? "identity" : "unresolved"),
            "symmetry_relation_is_meaningful" => resid <= 2.0,
        ))
    end
    sort!(rows; by = r -> -r["multiplicity_count"])
    return rows
end
