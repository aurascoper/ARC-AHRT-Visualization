# ==============================================================================
# src/structures/VesselDescriptors.jl
# Geometric descriptors computed from real coordinates. Every field name carries
# its unit; mmCIF Cartn_* are Angstrom, so lengths are reported in both A and nm.
#
# No descriptor is reported at a single threshold. The subunit contact graph is
# swept over a cutoff range and the whole curve is returned, because a single
# cutoff is how a threshold artifact gets mistaken for a structural fact.
# ==============================================================================

const ANGSTROM_PER_NM = 10.0

"""
    radius_of_gyration_angstrom(coords) -> Float64

Unweighted R_g over the supplied 3 x N coordinates (all atoms equal). Mass
weighting would need element masses; this is the geometric R_g and is labelled
as such.
"""
function radius_of_gyration_angstrom(coords::AbstractMatrix{<:Real})
    n = size(coords, 2)
    n == 0 && throw(ArgumentError("radius_of_gyration: no coordinates"))
    c = vec(sum(coords, dims=2)) ./ n
    acc = 0.0
    @inbounds for j in 1:n
        dx = coords[1, j] - c[1]; dy = coords[2, j] - c[2]; dz = coords[3, j] - c[3]
        acc += dx * dx + dy * dy + dz * dz
    end
    return sqrt(acc / n)
end

"""
    inertia_descriptors(coords) -> NamedTuple

Principal moments of the geometric gyration tensor and the standard asphericity
and acylindricity. Asphericity is 0 for a perfect sphere.
"""
function inertia_descriptors(coords::AbstractMatrix{<:Real})
    n = size(coords, 2)
    c = vec(sum(coords, dims=2)) ./ n
    X = coords .- c
    S = (X * X') ./ n                      # gyration tensor, A^2
    lam = sort(eigvals(Symmetric(S)); rev=true)
    l1, l2, l3 = lam[1], lam[2], lam[3]
    rg2 = l1 + l2 + l3
    asph = l1 - 0.5 * (l2 + l3)
    acyl = l2 - l3
    kappa2 = rg2 > 0 ? (asph^2 + 0.75 * acyl^2) / (rg2^2) : 0.0
    return (principal_moments_angstrom2 = [l1, l2, l3],
            asphericity_angstrom2 = asph,
            acylindricity_angstrom2 = acyl,
            relative_shape_anisotropy_dimensionless = kappa2)
end

"""
    radial_shell_profile(coords; nbins=120) -> Dict

Distance of every atom from the assembly centroid, binned. This is what settles
a claimed particle diameter against the deposited coordinates: `inner_radius`
is the 1st percentile of the radial distribution and `outer_radius` the 99th,
with the extremes also reported so a spike is visible rather than clipped.
"""
function radial_shell_profile(coords::AbstractMatrix{<:Real}; nbins::Int=120)
    n = size(coords, 2)
    c = vec(sum(coords, dims=2)) ./ n
    r = Vector{Float64}(undef, n)
    @inbounds for j in 1:n
        dx = coords[1, j] - c[1]; dy = coords[2, j] - c[2]; dz = coords[3, j] - c[3]
        r[j] = sqrt(dx * dx + dy * dy + dz * dz)
    end
    sort!(r)
    pct(p) = r[clamp(ceil(Int, p * n), 1, n)]
    rmin, rmax = r[1], r[end]
    edges = range(rmin, rmax; length=nbins + 1)
    counts = zeros(Int, nbins)
    @inbounds for v in r
        k = clamp(ceil(Int, (v - rmin) / max(rmax - rmin, eps()) * nbins), 1, nbins)
        counts[k] += 1
    end
    return Dict(
        "min_radius_angstrom" => rmin,
        "p01_radius_angstrom" => pct(0.01),
        "median_radius_angstrom" => pct(0.50),
        "p99_radius_angstrom" => pct(0.99),
        "max_radius_angstrom" => rmax,
        "inner_radius_nm" => pct(0.01) / ANGSTROM_PER_NM,
        "outer_radius_nm" => pct(0.99) / ANGSTROM_PER_NM,
        "max_radius_nm" => rmax / ANGSTROM_PER_NM,
        "outer_diameter_nm" => 2 * pct(0.99) / ANGSTROM_PER_NM,
        "max_diameter_nm" => 2 * rmax / ANGSTROM_PER_NM,
        "shell_thickness_nm" => (pct(0.99) - pct(0.01)) / ANGSTROM_PER_NM,
        "histogram_bin_edges_angstrom" => collect(edges),
        "histogram_counts" => counts,
        "n_atoms" => n,
        "units" => Dict("radius" => "angstrom and nm", "histogram_counts" => "atom count per bin"),
    )
end

# --- subunit contact graph ---------------------------------------------------

"""
    subunit_contact_graph(centroids, cutoff_angstrom) -> (nv, edges)

Nodes are polymer subunit centroids; an edge joins two whose centroids are
within `cutoff_angstrom`. Centroid distance is a coarse proxy for an interface
and is labelled as such; an all-atom minimum-distance criterion is the Release B
refinement.
"""
function subunit_contact_graph(centroids::Vector{Vector{Float64}}, cutoff_angstrom::Real)
    nv = length(centroids)
    edges = Tuple{Int,Int}[]
    c2 = cutoff_angstrom^2
    @inbounds for i in 1:nv-1, j in i+1:nv
        a = centroids[i]; b = centroids[j]
        dx = a[1]-b[1]; dy = a[2]-b[2]; dz = a[3]-b[3]
        (dx*dx + dy*dy + dz*dz) <= c2 && push!(edges, (i, j))
    end
    return nv, edges
end

"Connected components by union-find; returns (beta_0, component_id_per_vertex)."
function connected_components(nv::Int, edges::Vector{Tuple{Int,Int}})
    parent = collect(1:nv)
    find(x) = (while parent[x] != x; parent[x] = parent[parent[x]]; x = parent[x]; end; x)
    for (u, v) in edges
        ru, rv = find(u), find(v)
        ru != rv && (parent[ru] = rv)
    end
    roots = Dict{Int,Int}(); comp = zeros(Int, nv)
    for v in 1:nv
        r = find(v)
        comp[v] = get!(roots, r, length(roots) + 1)
    end
    return length(roots), comp
end

"beta_1 = E - V + beta_0 for a 1-complex (same convention as the existing repo code)."
function betti_numbers(nv::Int, edges::Vector{Tuple{Int,Int}})
    uniq = Set{Tuple{Int,Int}}()
    for (u, v) in edges
        push!(uniq, u < v ? (u, v) : (v, u))
    end
    b0, _ = connected_components(nv, collect(uniq))
    e = length(uniq)
    return b0, max(0, e - nv + b0), e
end

"""
    algebraic_connectivity(nv, edges) -> Float64

lambda_2 of the combinatorial graph Laplacian (Fiedler value). Zero when the
graph is disconnected. Dimensionless.
"""
function algebraic_connectivity(nv::Int, edges::Vector{Tuple{Int,Int}})
    nv <= 1 && return 0.0
    L = zeros(Float64, nv, nv)
    for (u, v) in edges
        L[u, v] -= 1.0; L[v, u] -= 1.0
        L[u, u] += 1.0; L[v, v] += 1.0
    end
    ev = sort(eigvals(Symmetric(L)))
    return length(ev) >= 2 ? max(0.0, ev[2]) : 0.0
end

"""
    contact_graph_sweep(centroids; cutoffs) -> Vector{Dict}

The whole beta_0 / beta_1 / lambda_2 curve against the contact cutoff. Reported
in full; no single cutoff is privileged.
"""
function contact_graph_sweep(centroids::Vector{Vector{Float64}};
                             cutoffs::AbstractVector{<:Real}=20.0:5.0:140.0)
    out = Dict{String,Any}[]
    for cut in cutoffs
        nv, edges = subunit_contact_graph(centroids, cut)
        b0, b1, e = betti_numbers(nv, edges)
        push!(out, Dict{String,Any}(
            "cutoff_angstrom" => float(cut),
            "cutoff_nm" => float(cut) / ANGSTROM_PER_NM,
            "n_vertices_count" => nv,
            "n_edges_count" => e,
            "mean_degree_dimensionless" => nv > 0 ? 2 * e / nv : 0.0,
            "beta_0_count" => b0,
            "beta_1_count" => b1,
            "euler_characteristic_count" => nv - e,
            "algebraic_connectivity_dimensionless" => algebraic_connectivity(nv, edges),
        ))
    end
    return out
end

# --- interface classes and virtual deletion ----------------------------------

"""
    classify_interfaces(centroids, cutoff_angstrom; ndigits=0) -> Dict{String,Vector{Int}}

Group contact-graph edges into classes by their centroid separation rounded to
`ndigits`. In a symmetric assembly the distinct separations correspond to the
distinct interface types (five-fold, three-fold, two-fold and so on), so this
recovers interface classes from geometry without hand annotation.
"""
function classify_interfaces(centroids::Vector{Vector{Float64}}, cutoff_angstrom::Real; ndigits::Int=0)
    nv, edges = subunit_contact_graph(centroids, cutoff_angstrom)
    classes = Dict{String,Vector{Int}}()
    for (k, (u, v)) in enumerate(edges)
        a = centroids[u]; b = centroids[v]
        d = sqrt((a[1]-b[1])^2 + (a[2]-b[2])^2 + (a[3]-b[3])^2)
        key = string(round(d; digits=ndigits))
        push!(get!(classes, key, Int[]), k)
    end
    return classes, nv, edges
end

"""
    interface_deletion_analysis(centroids, cutoff_angstrom; ndigits=0) -> Dict

The central Release A measurement. For each interface class, delete every edge
of that class and report the change in beta_0, beta_1, lambda_2 and the subunit
centroid R_g.

LIMITATION, carried in the output: this measures structural importance to GRAPH
CONNECTIVITY. It is not mechanical load bearing, which would require an
elastic-network or normal-mode analysis. Geometric integrity, graph connectivity
and mechanical stability are three distinct properties and are not merged here.
"""
function interface_deletion_analysis(centroids::Vector{Vector{Float64}}, cutoff_angstrom::Real;
                                     ndigits::Int=0)
    classes, nv, edges = classify_interfaces(centroids, cutoff_angstrom; ndigits=ndigits)
    b0_0, b1_0, e0 = betti_numbers(nv, edges)
    lam0 = algebraic_connectivity(nv, edges)
    coords0 = reduce(hcat, centroids)
    rg0 = radius_of_gyration_angstrom(coords0)

    rows = Dict{String,Any}[]
    for (key, idxs) in sort(collect(classes); by = x -> -length(x[2]))
        keep = [edges[k] for k in 1:length(edges) if !(k in Set(idxs))]
        b0, b1, e = betti_numbers(nv, keep)
        lam = algebraic_connectivity(nv, keep)
        push!(rows, Dict{String,Any}(
            "interface_class_separation_angstrom" => parse(Float64, key),
            "interface_class_separation_nm" => parse(Float64, key) / ANGSTROM_PER_NM,
            "n_edges_in_class_count" => length(idxs),
            "fraction_of_all_edges_dimensionless" => length(idxs) / max(1, e0),
            "delta_beta_0_count" => b0 - b0_0,
            "delta_beta_1_count" => b1 - b1_0,
            "delta_algebraic_connectivity_dimensionless" => lam - lam0,
            "beta_0_after_count" => b0,
            "beta_1_after_count" => b1,
            "algebraic_connectivity_after_dimensionless" => lam,
            "fragments_assembly" => b0 > b0_0,
        ))
    end

    return Dict(
        "cutoff_angstrom" => float(cutoff_angstrom),
        "cutoff_nm" => float(cutoff_angstrom) / ANGSTROM_PER_NM,
        "baseline" => Dict(
            "n_subunits_count" => nv,
            "n_edges_count" => e0,
            "beta_0_count" => b0_0,
            "beta_1_count" => b1_0,
            "algebraic_connectivity_dimensionless" => lam0,
            "subunit_centroid_radius_of_gyration_angstrom" => rg0,
            "subunit_centroid_radius_of_gyration_nm" => rg0 / ANGSTROM_PER_NM,
        ),
        "interface_classes" => rows,
        "n_interface_classes_count" => length(classes),
        "limitation" => "Measures structural importance to graph connectivity only. " *
                        "NOT mechanical load bearing: that requires an elastic-network or normal-mode analysis. " *
                        "Geometric integrity, graph connectivity and mechanical stability are three distinct properties.",
        "units" => Dict("separation" => "angstrom and nm",
                        "delta_beta_0_count" => "count",
                        "delta_beta_1_count" => "count",
                        "algebraic_connectivity" => "dimensionless (Laplacian eigenvalue)"),
    )
end

# --- top level ---------------------------------------------------------------

"""
    describe_assembly(asm; cutoffs) -> Dict

All descriptors for one generated assembly, with units on every field.
"""
function describe_assembly(asm::Assembly; cutoffs::AbstractVector{<:Real}=20.0:5.0:140.0)
    poly = [c for c in asm.chains if c.is_polymer]
    isempty(poly) && throw(ArgumentError("assembly $(asm.entry_id) has no polymer chains"))
    allcoords = reduce(hcat, (c.coords for c in poly))
    centroids = [c.centroid for c in poly]
    inert = inertia_descriptors(allcoords)
    rg = radius_of_gyration_angstrom(allcoords)

    return Dict(
        "entry_id" => asm.entry_id,
        "title" => asm.title,
        "source_sha256" => asm.source_sha256,
        "resolution_angstrom" => isnan(asm.resolution_angstrom) ? nothing : asm.resolution_angstrom,
        "point_symmetry_declared" => asm.point_symmetry_declared,
        "n_polymer_chains_count" => asm.n_polymer_chains,
        "n_chain_instances_count" => length(asm.chains),
        "declared_oligomeric_count" => asm.declared_oligomeric_count,
        "n_atoms_in_assembly_count" => size(allcoords, 2),
        "radius_of_gyration_angstrom" => rg,
        "radius_of_gyration_nm" => rg / ANGSTROM_PER_NM,
        "principal_moments_angstrom2" => inert.principal_moments_angstrom2,
        "asphericity_angstrom2" => inert.asphericity_angstrom2,
        "relative_shape_anisotropy_dimensionless" => inert.relative_shape_anisotropy_dimensionless,
        "radial_profile" => radial_shell_profile(allcoords),
        "contact_graph_sweep" => contact_graph_sweep(centroids; cutoffs=cutoffs),
        "triangulation" => derive_triangulation(asm.n_polymer_chains;
            declared_oligomeric_count = asm.declared_oligomeric_count),
    )
end
