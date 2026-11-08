# M7 / 6TAP Level-2 Symmetry-Expanded Residue Interface Atlas Specification

**Document ID**: `M7-ARC-LEVEL2-RESIDUE-SPEC-v1.0.0`
**Target Structure**: `PDB 6TAP` / `EMD-10423` (*Drosophila melanogaster* `dArc1` capsid biological assembly `1`, `3.50 A` cryo-EM reconstruction; Erlendsson et al., *Nat. Neurosci.* 2020, DOI: `10.1038/s41593-019-0569-y`)
**Schema Contract**: `schemas/m7_arc_level2_residue_interface_graph.schema.json`
**Reference Implementation**: `src/m7_arc_level2_interface_atlas.py`
**Validation Receipt**: `artifacts/m7_arc_level2_validation_receipt.json`
**Evidence Class**: `STRUCTURE_DERIVED_GEOMETRIC_CENSUS`

---

## 1. Primary Level-2 Object & Epistemic Separation from Level 1

### 1.1 Primary Level-2 Object: Symmetry-Expanded Residue-Interface Graph

The foundational mathematical and computational object of Level 2 is the **symmetry-expanded residue-interface graph**:

$$\mathcal{G}_{\text{L2}}(d_{\text{cut}}) = \left(\mathcal{V}_{\text{res}}, \mathcal{E}_{\text{contact}}(d_{\text{cut}})\right)$$

constructed directly from the deposited atomic coordinates (`_atom_site`) and the 60 icosahedral symmetry operators (`_pdbx_struct_oper_list`, `oper_expression = "(1-60)"`, icosahedral point group `I`, `_em_single_particle_entity.point_symmetry = "I"`, biological assembly `1`, RCSB global symmetry `I`; note that cubic space-group notation `I432` No. 211 is octahedral and is not applicable to `6TAP`) of `PDB 6TAP` biological assembly `1`.

Unlike subunit-centroid approximations or unsupervised feature clustering, $\mathcal{G}_{\text{L2}}(d_{\text{cut}})$ is a deterministic, symmetry-equivariant, mmCIF-native geometric graph in which:
- Every vertex $v \in \mathcal{V}_{\text{res}}$ is an explicitly indexed residue (or ligand/ion) symmetry instance with preserved mmCIF chemical entity, label/author chain, operator, and sequence identifiers.
- Every edge $e \in \mathcal{E}_{\text{contact}}(d_{\text{cut}})$ represents a verified inter-subunit heavy-atom contact at cutoff $d_{\text{cut}} \in \{3.5, 4.0, 4.5, 5.0, 5.5, 6.0\}\text{ \AA}$, annotated with exact minimum heavy-atom Euclidean distance, heavy-atom pair hit multiplicity, solvent-accessible surface area (SASA) burial, chemical interaction candidates, and exact $\text{SO}(3)$ operator-to-axis symmetry relations.

### 1.2 Explicit Separation from Level 1 Exploratory K-Means (`v1` and `v2`)

Level 1 (`src/m7_arc_cluster_bridge.py` `v1` and `src/m7_arc_cluster_bridge_v2.py` `v2`) applied unsupervised k-means clustering to per-subunit summary feature vectors over the 240 `CA` instances (`evidence_class = CALCULATED_EXPLORATORY`).

1. **Level 1 is NOT the foundation for Level 2**:
   - Level 1 `v1` (`k=3`, silhouette `0.248`, seed stability `0.649`, cluster sizes `91/85/64`) and Level 1 `v2` (`k=2`, silhouette `0.332`, cluster sizes `60/180`, label `CHAIN_ENVIRONMENT_CLASS_B_VS_ACD`) are retained strictly as historical exploratory coordinate-feature descriptors (`PRESERVED_AS_HISTORICAL_EXPLORATORY_EVIDENCE_NOT_FOUNDATION_FOR_LEVEL2`).
   - Level-2 residue annotations, interface classes, buried surface areas, and symmetry-axis assignments are derived **exclusively** from the mmCIF coordinates and `_pdbx_struct_oper_list` matrices of `6TAP.cif`, never conditioned on or filtered by Level-1 k-means cluster labels.
2. **Correction of Level-1 `v2` Follow-up Rationale on Homotypic Contacts (`:Superseded`)**:
   - The earlier Level-1 `v2` follow-up rationale claiming that Chain B separates from `A/C/D` because `B` alone lacks same-class contacts is marked `:Superseded`: **"B and D both lack homotypic interfaces in the 6TAP residue-contact census (B-B = 0 and D-D = 0 across 3.5..6.0 A; only A-A = 41 and C-C = 7 are homotypic at 5.0 A)"** (alongside `A-C = 0`). Only Chain A (`A-A` across the 5-fold `C5` pentamer ring, `op_id = 36/57`) and Chain C (`C-C` across the 3-fold `C3` hexamer ring, `op_id = 33/41`) form homotypic self-contacts.
   - **Level-2 Re-Evaluation of the `B`-versus-`A/C/D` `k=2` Partition**: In the Level-2 symmetry-expanded interface graph, `B` separates from `A/C/D` not by unique absence of homotypic contacts (since `D-D = 0` as well), but by its heterotypic partner-composition and interface-area profile: `B` forms two `B-C` interfaces (`1,272.29` and `1,365.72 A^2`) and two `B-D` interfaces (`1,342.37` and `344.66 A^2`) with only a marginal `48.65 A^2` `A-B` contact (`1` residue pair at `4.6059 A`, vanishing at `<= 4.5 A`), whereas `D` anchors the pentamer-to-hexamer seam via `A-D_op1` (`1,287.13 A^2`, `27` pairs) and `A-D_op36` (`179.98 A^2`, `5` pairs) alongside `C-D_op1` (`1,487.53 A^2`, `34` pairs).

---

## 2. Frozen Node & Edge Identity Contracts

### 2.1 Frozen Node Identity

Every node $v \in \mathcal{V}_{\text{res}}$ in the Level-2 symmetry-expanded residue graph is uniquely identified by the frozen 8-component tuple:

$$\text{NodeID}(v) = \texttt{structure\_id} \cdot \texttt{assembly\_id} \cdot \texttt{entity\_id} \cdot \texttt{asym\_id} \cdot \texttt{operator\_id} \cdot \texttt{label\_seq\_id} \cdot \texttt{auth\_seq\_id} \cdot \texttt{insertion\_code}$$

serialized in canonical dot-delimited form:
`{structure_id}.{assembly_id}.{entity_id}.{asym_id}.{operator_id}.{label_seq_id}.{auth_seq_id}.{insertion_code}`

| Component | Source mmCIF Field | Type | Value Domain in `6TAP` |
| :--- | :--- | :--- | :--- |
| `structure_id` | `_entry.id` | `String` | `"6TAP"` |
| `assembly_id` | `_pdbx_struct_assembly.id` | `String` | `"1"` |
| `entity_id` | `_atom_site.label_entity_id` | `String` | `"1"` (`dArc1` polypeptide strands `A..E`) or `"2"` (`ZN` non-polymer ligand strand `F`); paired with `functional_component_id` (`"CA_CAPSID_SHELL"`, `"NC_ZINC_KNUCKLE_FRAGMENT"`, `"ZN_COORDINATION_LIGAND"`) and `biological_domain_or_fragment_designation` |
| `asym_id` | `_atom_site.label_asym_id` | `String` | `"A"`, `"B"`, `"C"`, `"D"` (`CA` capsid shell); `"E"` (`NC` zinc-knuckle fragment of Entity `1`); `"F"` (`ZN` ligand of Entity `2`, where `_atom_site.auth_asym_id = "E"`) |
| `operator_id` | `_pdbx_struct_oper_list.id` | `String` | `"1"` through `"60"` |
| `label_seq_id` | `_atom_site.label_seq_id` | `Int64 \| null` | `42..205` (`A`), `41..204` (`B`), `41..205` (`C`), `41..206` (`D`), `224..250` (`E`), `null` / `"."` (`F`) |
| `auth_seq_id` | `_atom_site.auth_seq_id` | `Int64` | `42..205` (`A`), `41..204` (`B`), `41..205` (`C`), `41..206` (`D`), `224..250` (`E`), `301` (`F`, `_pdbx_nonpoly_scheme.auth_seq_num = 251`) |
| `insertion_code` | `_atom_site.pdbx_PDB_ins_code` | `String` | `"."` (no insertion codes in `6TAP`, preserved verbatim when present in general mmCIF inputs) |

### 2.2 Frozen Edge Identity

Every edge $e \in \mathcal{E}_{\text{contact}}(d_{\text{cut}})$ is uniquely identified by the frozen 4-component tuple:

$$\text{EdgeID}(e) = \texttt{residue\_instance\_i} \cdot \texttt{residue\_instance\_j} \cdot \texttt{cutoff\_definition\_id} \cdot \texttt{contact\_metric\_id}$$

serialized in canonical double-underscore-delimited form:
`{residue_instance_i}__{residue_instance_j}__{cutoff_definition_id}__{contact_metric_id}`

where:
- `residue_instance_i` and `residue_instance_j` are canonical 8-component `NodeID` strings of distinct subunit instances ($\text{instance}(i) \ne \text{instance}(j)$). For the asymmetric-unit generator representation, `residue_instance_i` resides in `operator_id = "1"` and `residue_instance_j` resides in `operator_id = op_j` with canonical chain ordering $\text{asym\_id}_i \le \text{asym\_id}_j$.
- `cutoff_definition_id` $\in \{\text{"HA\_3.5A"}, \text{"HA\_4.0A"}, \text{"HA\_4.5A"}, \text{"HA\_5.0A"}, \text{"HA\_5.5A"}, \text{"HA\_6.0A"}\}$.
- `contact_metric_id` $= \text{"MIN\_HEAVY\_ATOM\_EUCLIDEAN"}$.

### 2.3 Frozen Edge Attributes

Every edge in $\mathcal{E}_{\text{contact}}(d_{\text{cut}})$ carries the following 8 mandatory attributes:

1. `minimum_heavy_atom_distance` (`Float64`, $\text{\AA}$): Exact minimum Euclidean distance $\min_{a \in r_i, b \in r_j} \|\mathbf{x}_a - \mathbf{x}_b\|_2$ over all heavy atoms of the two residue instances.
2. `atom_pair_count` (`Int64`): Count of heavy-atom pairs $(a \in r_i, b \in r_j)$ satisfying $\|\mathbf{x}_a - \mathbf{x}_b\|_2 \le d_{\text{cut}}$.
3. `residue_pair_contact` (`Object`): Structured contact descriptor recording `is_contact = true`, `chain_pair`, `comp_i`, `seq_i`, `closest_atom_i`, `comp_j`, `seq_j`, `closest_atom_j`, and `symmetry_multiplicity_in_assembly = 60`.
4. `buried_surface_area` (`Object`, $\text{\AA}^2$): Structured SASA burial descriptor recording `direct_pair_occluded_sasa_A2` (surface area of $r_i$ directly occluded by atoms of $r_j$ plus surface area of $r_j$ directly occluded by atoms of $r_i$), `residue_i_interface_delta_sasa_A2`, and `residue_j_interface_delta_sasa_A2`.
5. `interaction_features` (`Object`): Structured chemical feature classification (`is_hydrophobic_contact_4p5A`, `is_hydrophobic_contact_5p0A`, `is_hbond_candidate_3p5A`, `is_salt_bridge_candidate_4p0A`, `is_pi_stacking_candidate_7p0A`, and `feature_tags`).
6. `interface_class` (`String`): Specific structural interface identifier combining chain pair, relative operator, and symmetry context (one of the 10 active capsid interface classes: `A-A_C5_pentamer_op36`, `A-B_C2_peripheral_op2`, `A-D_C1_intra_asu_op1`, `A-D_C5_inter_asu_seam_op36`, `B-C_C1_intra_asu_op1`, `B-C_C3_inter_asu_hexamer_op33`, `B-D_C2_dimer_op2`, `B-D_C3_inter_asu_vertex_op33`, `C-C_C3_trimer_op33`, `C-D_C1_intra_asu_op1`).
7. `operator_pair` (`Array[String, 2]`): Ordered operator pair `["1", op_j]` relating the two contacting subunit instances.
8. `symmetry_axis_relation` (`Object`): Exact operator-derived rotation relation containing `relative_operator_id`, `symmetry_symbol` (`"C1_identity"`, `"C2"`, `"C3"`, `"C5"`), `rotation_order` ($1, 2, 3, 5$), `rotation_angle_deg`, `axis_unit_vector_xyz`, and `is_true_crystallographic_self_rotation` (`true` iff $\text{asym\_id}_i = \text{asym\_id}_j$).

---

## 3. mmCIF Entity, Asym, and Operator Inventory Reconciliation

Direct parsing of `structures/6TAP.cif` (`SHA-256 = 61397d5502c44f6fd83d887de353a698ed6f70fbd00604e2fc1c85fd1cd68e69`, `649,219` bytes) establishes the exact reconciliation between deposited chemical entities (`_entity`, `_entity_poly`), asymmetric-unit strands (`_struct_asym`), atomic records (`_atom_site`), and biological assembly instances (`_pdbx_struct_assembly_gen`):

### 3.1 Deposited `_entity`, `_entity_poly`, `_struct_asym`, and `_atom_site.label_entity_id` Reconciliation for Chains `A..F`

In `6TAP.cif`, the `_entity` category defines **2 chemical entities**:
- `_entity.id = 1` (`polymer`, `polypeptide(L)`, `man`, `'Activity-regulated cytoskeleton associated protein 1'`, `dArc1`, `UniProt Q7K1U0`, `254` aa, formula weight `28,921.201 Da`, `_entity.pdbx_number_of_molecules = 5`, `_entity_poly.pdbx_strand_id = A,B,C,D,E`): Full-length 254-residue polypeptide sequence. Strands `A, B, C, D` resolve the 4 quasi-equivalent capsid (`CA`) domains (`residues 41..206`), whereas strand `E` (`residues 224..250`, with the flexible linker `207..223` unmodeled) is a resolved C-terminal nucleocapsid (`NC`) segment of the **same** deposited `dArc1` polypeptide Entity `1` (`label_entity_id = "1"`, `struct_asym_entity_id = "1"`), NOT an independently deposited polymer `"Entity 2"`.
- `_entity.id = 2` (`non-polymer`, `syn`, `'ZINC ION'`, `comp_id = ZN`, formula weight `65.409 Da`, `_entity.pdbx_number_of_molecules = 1`, `_pdbx_nonpoly_scheme.asym_id = F`, `_pdbx_nonpoly_scheme.pdb_strand_id = E`, `pdb_seq_num = 301`, `auth_seq_num = 251`, `_atom_site.label_seq_id = .`): The `ZN` ion is non-polymer ligand Entity `2` (`label_entity_id = "2"`, `struct_asym_entity_id = "2"`), NOT `"Entity 3"`.

To distinguish the `CA` capsid shell, the `NC` zinc-knuckle fragment, and the `ZN` coordination ligand without ever mislabeling `NC` as `"Entity 2"` or `ZN` as `"Entity 3"`, Level 2 records `functional_component_id` (`"CA_CAPSID_SHELL"`, `"NC_ZINC_KNUCKLE_FRAGMENT"`, `"ZN_COORDINATION_LIGAND"`) and `biological_domain_or_fragment_designation` alongside the literal mmCIF `label_entity_id` and `struct_asym_entity_id` (`"1"` for `A..E`, `"2"` for `F`):

| `label_asym_id` | `auth_asym_id` | `label_entity_id` | `struct_asym_entity_id` | `functional_component_id` | `biological_domain_or_fragment_designation` | Group | ASU Residues | ASU Atoms | `label_seq_id` Range | `auth_seq_id` Range | Mean $B_{\text{iso}}$ ($\text{\AA}^2$) | Operators | Assembly Instances |
| :---: | :---: | :---: | :---: | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `A` | `A` | `"1"` | `"1"` | `CA_CAPSID_SHELL` | `CA_capsid_shell_domain_pentamer_A_41_205` | `ATOM` | 164 | 1,343 | `42..205` | `42..205` | $76.63 \pm 7.68$ | 60 | 60 |
| `B` | `B` | `"1"` | `"1"` | `CA_CAPSID_SHELL` | `CA_capsid_shell_domain_hexamer_B_43_205` | `ATOM` | 164 | 1,339 | `41..204` | `41..204` | $109.90 \pm 8.81$ | 60 | 60 |
| `C` | `C` | `"1"` | `"1"` | `CA_CAPSID_SHELL` | `CA_capsid_shell_domain_hexamer_C_43_206` | `ATOM` | 165 | 1,347 | `41..205` | `41..205` | $112.53 \pm 10.61$ | 60 | 60 |
| `D` | `D` | `"1"` | `"1"` | `CA_CAPSID_SHELL` | `CA_capsid_shell_domain_hexamer_D_41_206` | `ATOM` | 166 | 1,356 | `41..206` | `41..206` | $112.17 \pm 10.17$ | 60 | 60 |
| `E` | `E` | `"1"` | `"1"` | `NC_ZINC_KNUCKLE_FRAGMENT` | `NC_zinc_knuckle_C_terminal_fragment_224_250_of_entity_1` | `ATOM` | 27 | 233 | `224..250` | `224..250` | $91.58 \pm 32.25$ | 60 | 60 |
| `F` | `E` | `"2"` | `"2"` | `ZN_COORDINATION_LIGAND` | `ZN_tetrahedral_CCHC_ligand_entity_2` | `HETATM` | 1 | 1 | `null` (`.`) | `301..301` | $93.87 \pm 0.00$ | 60 | 60 |
| **Total** | -- | -- | -- | **All Deposited Asym Units** | -- | -- | **687** | **5,619** | -- | -- | -- | **60** | **360** |

**Inventory Reconciliation Summary**:
- `4` deposited `CA` chains (`A, B, C, D`, `659` residues, `5,385` heavy atoms) $\times 60$ operators $= 240$ `CA` instances (`39,540` `CA` residue nodes, `323,100` `CA` atoms; Caspar-Klug $T = 240 / 60 = 4$).
- `1` deposited `NC` peptide chain (`E`, `27` residues, `233` heavy atoms) $\times 60$ operators $= 60$ `NC` symmetry copies (`1,620` `NC` residue nodes, `13,980` `NC` atoms).
- `240` `CA` instances $+ 60$ `NC` instances $= 300$ polymer chain instances, reconciling `_pdbx_struct_assembly.oligomeric_count = 300` (`300-meric`).
- `1` deposited `Zn` ion (`F`, `1` atom) $\times 60$ operators $= 60$ `Zn` symmetry copies (`360` total chain/asym instances, `41,220` total residue/ion nodes, `337,140` total atoms).
- `0` water (`HOH`) molecules and `0` other ligands are present in `6TAP.cif`.

---

## 4. Work Package Specifications & Exact Verified Tables

### 4.1 Work Package `L2-A`: Multi-Cutoff Residue-Pair Contact Census (`3.5 A` to `6.0 A`)

#### 4.1.1 Three Explicit Contact-Count Definitions

To eliminate ambiguity between asymmetric-unit residue pairs, assembly-wide residue-instance edges, and atomic pair hits, Level 2 defines and reports three distinct integer metrics for every chain-class pair $C_1\text{-}C_2$:

1. **`unique_asym_residue_pair_count`**:
   The cardinality of the set of distinct unordered asymmetric-unit residue pairs $\{(c_1, \text{seq}_1), (c_2, \text{seq}_2)\}$ with $c_1 \in C_1, c_2 \in C_2$ that exhibit at least one cross-instance heavy-atom pair within distance $d_{\text{cut}}$ under at least one relative symmetry operator $\text{op}_k \in \{1..60\}$.
2. **`instance_residue_pair_observation_count`**:
   The cardinality of the set of distinct unordered residue-instance edges $\{((\text{op}_i, c_1, \text{seq}_1), (\text{op}_j, c_2, \text{seq}_2))\}$ with $(\text{op}_i, c_1) \ne (\text{op}_j, c_2)$ across all $240$ `CA` subunit instances in the full icosahedral assembly. Because the 60 operators of $I$ act freely and transitively on the 60 asymmetric units without fixed points on any subunit centroid, every unique asymmetric-unit interface contact has exact symmetry multiplicity $60$:
   $$\texttt{instance\_residue\_pair\_observation\_count} = 60 \times \texttt{unique\_asym\_residue\_pair\_count}$$
3. **`heavy_atom_pair_hit_count`**:
   The total number of distinct unordered cross-instance heavy-atom pairs $(a, b)$ satisfying $\|\mathbf{x}_a - \mathbf{x}_b\|_2 \le d_{\text{cut}}$ across the full 240-subunit assembly. Multiple atom pairs between the same two residues increment `heavy_atom_pair_hit_count` while representing a single residue-residue contact edge.

#### 4.1.2 Complete Multi-Cutoff Census Across All 10 Chain-Class Pairs (`3.5 A` to `6.0 A`)

**Metric 1: `unique_asym_residue_pair_count`**

| Cutoff ($\text{\AA}$) | `A-A` | `A-B` | `A-C` | `A-D` | `B-B` | `B-C` | `B-D` | `C-C` | `C-D` | `D-D` | Total Unique Asym Pairs |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **3.5** | 11 | 0 | 0 | 8 | 0 | 16 | 5 | 2 | 8 | 0 | **50** |
| **4.0** | 21 | 0 | 0 | 12 | 0 | 24 | 8 | 2 | 12 | 0 | **79** |
| **4.5** | 29 | 0 | 0 | 21 | 0 | 37 | 23 | 4 | 22 | 0 | **136** |
| **5.0** | **41** | **1** | **0** | **32** | **0** | **54** | **35** | **7** | **34** | **0** | **204** |
| **5.5** | 55 | 1 | 0 | 40 | 0 | 77 | 53 | 10 | 48 | 0 | **284** |
| **6.0** | 72 | 4 | 0 | 52 | 0 | 108 | 67 | 11 | 57 | 0 | **371** |

**Metric 2: `instance_residue_pair_observation_count` ($= 60 \times \texttt{unique\_asym\_residue\_pair\_count}$)**

| Cutoff ($\text{\AA}$) | `A-A` | `A-B` | `A-C` | `A-D` | `B-B` | `B-C` | `B-D` | `C-C` | `C-D` | `D-D` | Total Instance-Residue Edges |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **3.5** | 660 | 0 | 0 | 480 | 0 | 960 | 300 | 120 | 480 | 0 | **3,000** |
| **4.0** | 1,260 | 0 | 0 | 720 | 0 | 1,440 | 480 | 120 | 720 | 0 | **4,740** |
| **4.5** | 1,740 | 0 | 0 | 1,260 | 0 | 2,220 | 1,380 | 240 | 1,320 | 0 | **8,160** |
| **5.0** | **2,460** | **60** | **0** | **1,920** | **0** | **3,240** | **2,100** | **420** | **2,040** | **0** | **12,240** |
| **5.5** | 3,300 | 60 | 0 | 2,400 | 0 | 4,620 | 3,180 | 600 | 2,880 | 0 | **17,040** |
| **6.0** | 4,320 | 240 | 0 | 3,120 | 0 | 6,480 | 4,020 | 660 | 3,420 | 0 | **22,260** |

**Metric 3: `heavy_atom_pair_hit_count` (Symmetry-Equivariant 240-Subunit Assembly Census)**

| Cutoff ($\text{\AA}$) | `A-A` | `A-B` | `A-C` | `A-D` | `B-B` | `B-C` | `B-D` | `C-C` | `C-D` | `D-D` | Total Heavy-Atom Pair Hits |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **3.5** | 780 | 0 | 0 | 1,200 | 0 | 3,360 | 360 | 120 | 840 | 0 | **6,660** |
| **4.0** | 2,400 | 0 | 0 | 3,120 | 0 | 6,300 | 900 | 120 | 1,500 | 0 | **14,340** |
| **4.5** | 5,340 | 0 | 0 | 6,180 | 0 | 12,660 | 3,240 | 420 | 4,020 | 0 | **31,860** |
| **5.0** | **13,080** | **180** | **0** | **12,660** | **0** | **22,380** | **7,980** | **1,530** | **8,580** | **0** | **66,390** |
| **5.5** | 24,480 | 480 | 0 | 21,480 | 0 | 38,700 | 16,920 | 3,390 | 15,480 | 0 | **120,930** |
| **6.0** | 41,340 | 1,260 | 0 | 34,860 | 0 | 62,160 | 29,760 | 5,880 | 26,760 | 0 | **202,020** |

*Note on 6-decimal operator rounding at boundary cutoffs*: Because `_pdbx_struct_oper_list` matrices in `6TAP.cif` are deposited to 6 decimal places ($\sim 10^{-6}$ precision), a heavy-atom pair sitting within $\sim 10^{-4}\text{ \AA}$ of an integer or half-integer cutoff (such as `C:GLN197:NE2 -- C:HIS186:CE1` at $d = 4.999979\text{ \AA}$ under `op_33` vs $5.000098\text{ \AA}$ under its inverse `op_41`) can straddle the cutoff across symmetry copies. The table above reports the canonical symmetry-equivariant count generated from reference `operator_id = "1"` ($\times 60$ for heterotypic pairs; $\times 60 / 2$ over directed `op_1 -> op_k` pairs for homotypic `A-A` and `C-C`). Evaluating all $\binom{240}{2}$ transformed copies independently without reference-operator canonicalization yields `31,832` (`4.5 A`), `66,337` (`5.0 A`), `120,938` (`5.5 A`), and `202,016` (`6.0 A`), while Metrics 1 and 2 (`204` and `12,240` at `5.0 A`) are invariant under both conventions.

#### 4.1.3 Decomposition into the 10 Symmetry-Inequivalent Sub-Interfaces at `5.0 A`

Because heterotypic pairs `A-D`, `B-C`, and `B-D` each span **two** distinct relative symmetry operators (`op_id = 1` intra-ASU vs `op_id = 36` or `op_id = 33` inter-ASU, plus `op_id = 2` for `B-D`), the 204 unique residue-pair contacts at `5.0 A` decompose into **10 symmetry-inequivalent sub-interfaces** (generating $10 \times 60 = 600$ undirected subunit-subunit edges across the 240-subunit assembly, i.e., $5.00$ partners per subunit):

| Interface Class | Chain Pair | Operator Pair | Inverse Op | Symmetry | Angle ($\text{deg}$) | Centroid Dist ($\text{\AA}$) | Min Atom Dist ($\text{\AA}$) | Unique Asym Pairs | Instance Edges | Heavy-Atom Hits |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `A-A_C5_pentamer_op36` | `A-A` | `["1", "36"]` | `"57"` | `C5` | 72.0 | 31.980 | 3.0278 | 41 | 2,460 | 13,080 |
| `A-B_C2_peripheral_op2` | `A-B` | `["1", "2"]` | `"2"` | `C2` | 180.0 | 64.326 | 4.6059 | 1 | 60 | 180 |
| `A-D_C1_intra_asu_op1` | `A-D` | `["1", "1"]` | `"1"` | `C1_identity` | 0.0 | 38.354 | 2.6868 | 27 | 1,620 | 11,100 |
| `A-D_C5_inter_asu_seam_op36` | `A-D` | `["1", "36"]` | `"57"` | `C5` | 72.0 | 63.911 | 1.9910 | 5 | 300 | 1,560 |
| `B-C_C1_intra_asu_op1` | `B-C` | `["1", "1"]` | `"1"` | `C1_identity` | 0.0 | 32.424 | 2.7119 | 26 | 1,560 | 7,260 |
| `B-C_C3_inter_asu_hexamer_op33` | `B-C` | `["1", "33"]` | `"41"` | `C3` | 120.0 | 40.004 | 1.8013 | 28 | 1,680 | 15,120 |
| `B-D_C2_dimer_op2` | `B-D` | `["1", "2"]` | `"2"` | `C2` | 180.0 | 34.429 | 3.0001 | 28 | 1,680 | 6,360 |
| `B-D_C3_inter_asu_vertex_op33` | `B-D` | `["1", "33"]` | `"41"` | `C3` | 120.0 | 65.085 | 4.0362 | 7 | 420 | 1,620 |
| `C-C_C3_trimer_op33` | `C-C` | `["1", "33"]` | `"41"` | `C3` | 120.0 | 65.504 | 3.2109 | 7 | 420 | 1,530 |
| `C-D_C1_intra_asu_op1` | `C-D` | `["1", "1"]` | `"1"` | `C1_identity` | 0.0 | 34.155 | 2.2642 | 34 | 2,040 | 8,580 |
| **Total (All 10 Sub-Interfaces)** | -- | -- | -- | -- | -- | -- | -- | **204** | **12,240** | **66,390** |

---

### 4.2 Work Package `L2-B`: Solvent-Accessible Surface Area (SASA) & Buried Surface Area (BSA)

Per-residue and per-interface solvent-accessible surface areas are computed via the deterministic Lee-Richards / Shrake-Rupley rolling-probe quadrature (`probe_radius = 1.40 A`, `100` Fibonacci golden-spiral points per atom, Bondi van der Waals radii: $r_{\text{C}} = 1.70\text{ \AA}$, $r_{\text{N}} = 1.55\text{ \AA}$, $r_{\text{O}} = 1.52\text{ \AA}$, $r_{\text{S}} = 1.80\text{ \AA}$, $r_{\text{Zn}} = 1.39\text{ \AA}$).

#### 4.2.1 Subunit-Level SASA Summary (`operator_id = "1"` Reference vs Full 240-Subunit Assembly)

| Chain | Isolated Subunit SASA ($\text{\AA}^2$) | Assembly SASA ($\text{\AA}^2$) | Total Subunit $\Delta\text{SASA}$ ($\text{\AA}^2$) | Buried Surface Fraction |
| :---: | :---: | :---: | :---: | :---: |
| `A` | 12,303.04 | 9,686.54 | 2,616.49 | 0.2127 (21.27%) |
| `B` | 12,362.50 | 10,203.19 | 2,159.32 | 0.1747 (17.47%) |
| `C` | 12,544.53 | 10,195.74 | 2,348.78 | 0.1872 (18.72%) |
| `D` | 12,390.95 | 10,066.97 | 2,323.98 | 0.1876 (18.76%) |

#### 4.2.2 Pairwise Buried Surface Area by Sub-Interface, Chain-Class Pair, and Whole-Capsid Multiplicity

For each contacting subunit pair $(U, V)$, $\text{BSA}_{\text{total}}(U, V) = \text{SASA}(U) + \text{SASA}(V) - \text{SASA}(U \cup V)$ and $\text{BSA}_{\text{per\_side}}(U, V) = \frac{1}{2}\text{BSA}_{\text{total}}(U, V)$. To prevent conflating the sum over the 10 representative symmetry-inequivalent sub-interfaces with the whole-capsid buried surface area, both the schema and receipt explicitly record:
- `representative_subinterface_pairwise_bsa_sum_A2 = 9570.22` (sum over the 10 representative symmetry-inequivalent sub-interfaces per asymmetric unit, NOT the whole-capsid buried area)
- `representative_subinterface_per_side_bsa_sum_A2 = 4785.10`
- `unordered_capsid_edge_multiplicity` per sub-interface (`60` for non-involution or heterotypic edges; `30` for homotypic `C2` involutions if any; in `6TAP`, all 10 active sub-interfaces are either heterotypic or homotypic `C5`/`C3` non-involutions, so every sub-interface has `unordered_capsid_edge_multiplicity = 60`)
- `whole_capsid_total_buried_surface_area_A2 = 574213.20` (`574,213.20 A^2` two-sided total across all `60` asymmetric units / `600` unordered subunit-subunit edges) and `whole_capsid_one_side_interface_area_A2 = 287106.60` (`287,106.60 A^2 = 60 * 4,785.10 A^2` one-side interface area across the whole capsid) alongside per-subinterface `whole_capsid_one_side_bsa_A2` (`60 * per_side_bsa_A2`) so `9,570.22 A^2` can never be mistaken for the whole-capsid buried area.

| Interface Class | Chain Pair | Operator Pair | `unordered_capsid_edge_multiplicity` | Subunit $i$ Buried ($\text{\AA}^2$) | Subunit $j$ Buried ($\text{\AA}^2$) | Total Pairwise BSA ($\text{\AA}^2$) | Per-Side BSA ($\text{\AA}^2$) | Whole-Capsid Two-Sided Sub-Interface BSA ($\text{\AA}^2$) | Whole-Capsid One-Side Sub-Interface BSA ($\text{\AA}^2$) | Chain-Pair Combined Total BSA ($\text{\AA}^2$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `A-A_C5_pentamer_op36` | `A-A` | `["1", "36"]` | 60 | 920.85 | 961.15 | **1,882.00** | 941.00 | 112,920.00 | 56,460.00 | **1,882.00** |
| `A-B_C2_peripheral_op2` | `A-B` | `["1", "2"]` | 60 | 34.07 | 14.58 | **48.65** | 24.33 | 2,919.00 | 1,459.80 | **48.65** |
| `A-D_C1_intra_asu_op1` | `A-D` | `["1", "1"]` | 60 | 645.07 | 642.06 | **1,287.13** | 643.56 | 77,227.80 | 38,613.60 | **1,467.11** (sum of `op1` + `op36`) |
| `A-D_C5_inter_asu_seam_op36` | `A-D` | `["1", "36"]` | 60 | 86.41 | 93.56 | **179.98** | 89.99 | 10,798.80 | 5,399.40 | (included in `A-D` `1,467.11`) |
| `B-C_C1_intra_asu_op1` | `B-C` | `["1", "1"]` | 60 | 650.07 | 622.22 | **1,272.29** | 636.15 | 76,337.40 | 38,169.00 | **2,638.01** (sum of `op1` + `op33`) |
| `B-C_C3_inter_asu_hexamer_op33` | `B-C` | `["1", "33"]` | 60 | 660.50 | 705.23 | **1,365.72** | 682.86 | 81,943.20 | 40,971.60 | (included in `B-C` `2,638.01`) |
| `B-D_C2_dimer_op2` | `B-D` | `["1", "2"]` | 60 | 666.69 | 675.68 | **1,342.37** | 671.18 | 80,542.20 | 40,270.80 | **1,687.03** (sum of `op2` + `op33`) |
| `B-D_C3_inter_asu_vertex_op33` | `B-D` | `["1", "33"]` | 60 | 145.73 | 198.93 | **344.66** | 172.33 | 20,679.60 | 10,339.80 | (included in `B-D` `1,687.03`) |
| `C-C_C3_trimer_op33` | `C-C` | `["1", "33"]` | 60 | 169.94 | 189.95 | **359.89** | 179.94 | 21,593.40 | 10,796.40 | **359.89** |
| `C-D_C1_intra_asu_op1` | `C-D` | `["1", "1"]` | 60 | 738.30 | 749.23 | **1,487.53** | 743.76 | 89,251.80 | 44,625.60 | **1,487.53** |
| `A-C`, `B-B`, `D-D` | -- | -- | 0 | 0.00 | 0.00 | **0.00** | 0.00 | 0.00 | 0.00 | **0.00** |
| **Representative Sum / Whole-Capsid Total** | -- | -- | **600 edges** | **4,717.63** | **4,852.59** | **`9,570.22`** (`representative_subinterface_pairwise_bsa_sum_A2`) | **`4,785.10`** (`representative_subinterface_per_side_bsa_sum_A2`) | **`574,213.20`** (`whole_capsid_total_buried_surface_area_A2`) | **`287,106.00` / `287,106.60`** (`whole_capsid_one_side_interface_area_A2`) | **`9,570.22`** |

---

### 4.3 Work Package `L2-C`: Candidate Interface Chemistry Classification

Because `PDB 6TAP` is a `3.50 A` cryo-EM heavy-atom coordinates model lacking experimental hydrogen positions, all `L2-C` non-covalent interaction counts are explicitly qualified as **geometric candidate counts** (`hydrophobic_candidate_pairs`, `pi_stacking_candidate_pairs`, `hbond_heavy_atom_distance_candidate_pairs`, `salt_bridge_heavy_atom_distance_candidate_pairs`) with mandatory geometric criteria and caveats:
- **`hydrophobic_candidate_pairs`**: Nonpolar side-chain carbon/sulfur atom pairs (`ALA, VAL, LEU, ILE, MET, PHE, TRP, PRO, TYR`, excluding backbone `N, CA, C, O, OXT`) within $d_{\text{C/S}} \le 4.5\text{ \AA}$ (and extended $d_{\text{C/S}} \le 5.0\text{ \AA}$).
- **`pi_stacking_candidate_pairs`**: Aromatic ring pairs (`PHE, TYR, TRP, HIS`) across active sub-interfaces with ring centroid distance `<= 7.0 A` aromatic ring-centroid distance with reported `ring_plane_angle_deg` (and alias `ring_normal_angle_deg` $\gamma_{\text{deg}} \in [0^\circ, 90^\circ]$) and `min_heavy_atom_distance_A` (8 of 10 candidates also satisfy minimum heavy-atom distance $\le 5.0\text{ \AA}$ and appear directly in the `5.0 A` edge list, while the two 3-fold vertex/trimer `HIS182-HIS182` candidates have $d_{\text{centroid}} = 6.808\text{ \AA}$ and $6.893\text{ \AA}$ with minimum heavy-atom distances $5.6214\text{ \AA}$ and $5.4825\text{ \AA}$).
- **`hbond_heavy_atom_distance_candidate_pairs`**: Donor-acceptor heavy-atom (`N, O, S`) pairs within `<= 3.5 A` N/O heavy-atom distance without explicit H-atom donor-H-acceptor angle (because hydrogens are not resolved at `3.50 A` cryo-EM resolution).
- **`salt_bridge_heavy_atom_distance_candidate_pairs`**: Cationic side-chain nitrogens (`ARG: NE, NH1, NH2`; `LYS: NZ`; `HIS: ND1, NE2`) to anionic carboxylate oxygens (`ASP: OD1, OD2`; `GLU: OE1, OE2`) within $d_{\text{NO}} \le 4.0\text{ \AA}$ heavy-atom distance without protonation-state or pKa assignment.

#### 4.3.1 Candidate Interface Chemistry Counts by Sub-Interface

| Interface Class | Chain Pair | Operator Pair | `hydrophobic_candidate_pairs` ($\le 4.5\text{ \AA}$) | Hydrophobic Extended ($\le 5.0\text{ \AA}$) | `pi_stacking_candidate_pairs` ($\le 7.0\text{ \AA}$) | `hbond_heavy_atom_distance_candidate_pairs` ($\le 3.5\text{ \AA}$) | `salt_bridge_heavy_atom_distance_candidate_pairs` ($\le 4.0\text{ \AA}$) | Key Conserved Motifs |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| `A-A_C5_pentamer_op36` | `A-A` | `["1", "36"]` | 2 | 4 | 0 | 4 | 3 | `ASP144-ARG56`, `ARG199-GLU117`, `TYR71-SER46` |
| `A-B_C2_peripheral_op2` | `A-B` | `["1", "2"]` | 0 | 0 | 0 | 0 | 0 | Single vdW grazing contact `LYS181-ASN204` (`4.61 A`) |
| `A-D_C1_intra_asu_op1` | `A-D` | `["1", "1"]` | 7 | 8 | 4 | 1 | 0 | C-lobe dimer: `PHE133-PHE133` (`4.22 A`, `35.3 deg`), `TYR126-PHE172` |
| `A-D_C5_inter_asu_seam_op36` | `A-D` | `["1", "36"]` | 0 | 0 | 0 | 2 | 0 | `ASN204(ND2)-HIS186(ND1/NE2)` (`2.66/2.91 A`) |
| `B-C_C1_intra_asu_op1` | `B-C` | `["1", "1"]` | 1 | 2 | 0 | 3 | 2 | `ARG56-ASP144` (`2.71/3.10 A`), `GLN127-HIS203` |
| `B-C_C3_inter_asu_hexamer_op33` | `B-C` | `["1", "33"]` | 6 | 8 | 4 | 5 | 0 | C-lobe dimer: `PHE133-PHE133` (`4.93 A`, `22.2 deg`), `ARG161-THR166` |
| `B-D_C2_dimer_op2` | `B-D` | `["1", "2"]` | 3 | 4 | 0 | 0 | 1 | `ASP144-ARG56` (`3.93 A`), N-lobe helix bundle |
| `B-D_C3_inter_asu_vertex_op33` | `B-D` | `["1", "33"]` | 0 | 0 | 1 | 0 | 0 | `HIS182-HIS182` (`6.81 A`, `72.4 deg`) |
| `C-C_C3_trimer_op33` | `C-C` | `["1", "33"]` | 0 | 0 | 1 | 1 | 0 | `HIS182-HIS182` (`6.89 A`, `85.5 deg`), `ASN205-LYS181` |
| `C-D_C1_intra_asu_op1` | `C-D` | `["1", "1"]` | 1 | 1 | 0 | 3 | 3 | `ARG56-ASP144` (`2.26/3.00 A`), `LYS84-ASP73` (`3.51 A`) |
| **Total** | -- | -- | **20** | **27** | **10** | **19** | **9** | -- |

---

### 4.4 Work Package `L2-D`: Exact Icosahedral Operator-to-Axis Mathematical Derivation

#### 4.4.1 Eigenvalue / Trace & Skew-Symmetric Axis Extraction

Instead of heuristic k-means clustering on subunit centroids (which introduced non-zero variance across exact symmetry copies in `m7_arc_level2_first.py`), `L2-D` derives the rotation angle $\theta_k$, rotation order $n_k$, and unit rotation axis $\mathbf{u}_k \in S^2$ directly from each operator matrix $R_k \in \text{SO}(3)$ ($k \in \{1..60\}$) in `_pdbx_struct_oper_list`:

1. **SO(3) Verification & Group Closure**:
   - Every deposited $3 \times 3$ matrix $R_k$ satisfies $\|R_k^\top R_k - I_{\text{3}}\|_{\max} \le 1.15 \times 10^{-6}$ and $|\det(R_k) - 1| \le 1.33 \times 10^{-6}$.
   - All $60 \times 60 = 3,600$ pairwise products $R_i R_j$ close onto a unique group element $R_m$ within max matrix-element residual $1.45 \times 10^{-6}$.
   - Least-squares solution of $(I - R_k)\mathbf{x}_0 = \mathbf{t}_k$ across $k = 2..60$ yields the unique global icosahedral fixed point (capsid center):
     $$\mathbf{x}_0 = [310.015983,\ 310.015978,\ 310.015954]^\top\text{ \AA}\quad (\text{max fixed-point residual } 3.325 \times 10^{-4}\text{ \AA})$$
2. **Rotation Angle & Order**:
   $$\theta_k = \arccos\!\left(\text{clip}\!\left(\frac{\text{tr}(R_k) - 1}{2}, -1, 1\right)\right) \in [0, \pi]$$
   - $\theta_k = 0^\circ$ ($\text{tr}(R_k) = 3$): order $n_k = 1$ (`C1_identity`, $1$ operator: `op_id = 1`).
   - $\theta_k \in \{72^\circ, 144^\circ\}$ ($\text{tr}(R_k) \in \{\frac{\sqrt{5}-1}{2}+1 \approx 1.618034,\ \frac{-\sqrt{5}-1}{2}+1 \approx -0.618034\}$): order $n_k = 5$ (`C5`, $24$ operators: $12$ with $\theta = 72^\circ$ and $12$ with $\theta = 144^\circ$, spanning **6 undirected 5-fold axis lines** / **12 directed vertex rays**).
   - $\theta_k = 120^\circ$ ($\text{tr}(R_k) = 0$): order $n_k = 3$ (`C3`, $20$ operators, spanning **10 undirected 3-fold axis lines** / **20 directed face-center rays**).
   - $\theta_k = 180^\circ$ ($\text{tr}(R_k) = -1$): order $n_k = 2$ (`C2`, $15$ involutory operators $R_k^2 = I$, spanning **15 undirected 2-fold axis lines** / **30 directed edge-midpoint rays**).
3. **Unit Rotation Axis $\mathbf{u}_k$**:
   - For $0^\circ < \theta_k < 180^\circ$ (`C5` and `C3`), the skew-symmetric part $R_k - R_k^\top = 2\sin\theta_k [\mathbf{u}_k]_\times$ directly determines the oriented unit axis:
     $$\mathbf{u}_k = \frac{1}{2\sin\theta_k}\begin{bmatrix} R_{k,32} - R_{k,23} \\ R_{k,13} - R_{k,31} \\ R_{k,21} - R_{k,12} \end{bmatrix}$$
   - For $\theta_k = 180^\circ$ (`C2`), $R_k$ is symmetric ($R_k = R_k^\top = 2\mathbf{u}_k\mathbf{u}_k^\top - I$), and $\mathbf{u}_k$ is the unit eigenvector corresponding to eigenvalue $\lambda = +1$ (canonicalized with the first non-zero coordinate positive).

#### 4.4.2 Active Contacting Operators in `6TAP`

| `operator_id` | Symmetry Symbol | Rotation Order | $\text{tr}(R_k)$ | $\theta_k$ ($\text{deg}$) | Unit Rotation Axis $\mathbf{u}_k = [u_x, u_y, u_z]$ | Translation Vector $\mathbf{t}_k$ ($\text{\AA}$) | Active Interfaces Generated from `op_id = 1` |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :--- |
| `"1"` | `C1_identity` | 1 | `3.000000` | `0.00` | `[0.000000, 0.000000, 1.000000]` | `[0.000, 0.000, 0.000]` | Intra-ASU quasi-equivalent: `A-D[op1]`, `B-C[op1]`, `C-D[op1]` |
| `"2"` | `C2` | 2 | `-1.000000` | `180.00` | `[0.425325, -0.309017, 0.850651]` | `[365.030, 805.302, 110.028]` | 2-fold involutory: `B-D[op2]`, `A-B[op2]`, `E-E[op2]` |
| `"33"` | `C3` | 3 | `0.000000` | `120.00` | `[-0.187593, 0.577350, -0.794654]` | `[61.345, 736.716, 520.774]` | 3-fold / quasi-6-fold ring: `B-C[op33]`, `B-D[op33]`, `C-C[op33]` |
| `"36"` | `C5` | 5 | `1.618034` | `72.00` | `[0.000000, 0.000000, 1.000000]` | `[509.059, -80.627, 0.000]` | 5-fold pentamer ring & seam: `A-A[op36]`, `A-D[op36]` |
| `"41"` | `C3` | 3 | `0.000000` | `120.00` | `[0.187593, -0.577350, 0.794654]` | `[798.061, 410.746, 110.028]` | Inverse of `op_33`: `C-B[op41]`, `D-B[op41]`, `C-C[op41]` |
| `"57"` | `C5` | 5 | `1.618034` | `72.00` | `[0.000000, 0.000000, -1.000000]` | `[-80.627, 509.059, 0.000]` | Inverse of `op_36`: `A-A[op57]`, `D-A[op57]` |

#### 4.4.3 Exact Subunit Radial Distance & Angular Alignment to Canonical Icosahedral Axes

Because all 60 operators are exact $\text{SO}(3)$ symmetries about $\mathbf{x}_0$, every symmetry copy of a given chain has **identical** angular distance ($\text{std} = 0.0000^\circ$) to the nearest operator-derived `C5`, `C3`, and `C2` axes:

| Chain (`asym_id`) | Domain Role | Radial Distance from $\mathbf{x}_0$ ($\text{\AA}$) | Nearest `C5` (5-Fold Vertex) Angle ($\text{deg}$) | Nearest `C3` (3-Fold / Quasi-6-Fold) Angle ($\text{deg}$) | Nearest `C2` (2-Fold Edge) Angle ($\text{deg}$) | Std Across 60 Copies ($\text{deg}$) |
| :---: | :--- | :---: | :---: | :---: | :---: | :---: |
| `A` | `CA` 5-fold pentamer subunit | 180.3827 | **8.6738** | 30.3810 | 23.0639 | 0.0000 |
| `B` | `CA` quasi-6-fold hexamer subunit | 173.0040 | 28.2676 | **11.4279** | 11.9972 | 0.0000 |
| `C` | `CA` quasi-6-fold hexamer subunit | 171.1227 | 27.7403 | 12.7680 | **10.7429** | 0.0000 |
| `D` | `CA` quasi-6-fold hexamer subunit | 173.6599 | 20.8574 | 22.8203 | **10.8940** | 0.0000 |
| `E` | `NC` lumenal zinc-knuckle peptide | 163.6982 | 29.0144 | 20.8869 | **2.7100** | 0.0000 |
| `F` | `Zn` coordination ion | 164.9754 | 29.2751 | 20.2633 | **2.5742** | 0.0000 |

#### 4.4.4 Proof that `operator_id` (`1..60`) is Deposit-Specific

In mmCIF files, `_pdbx_struct_oper_list.id` is an arbitrary row key assigned at deposition time, NOT a universal icosahedral group element label:
- In `6TAP.cif`, `op_id = "2"` is a 2-fold (`C2`, $180^\circ$) operator, `op_id = "3"` is a 3-fold (`C3`, $120^\circ$) operator, and `op_id = "36"` is the primary 5-fold (`C5`, $72^\circ$) pentamer generator.
- In `1OHF.cif` (Nudaurelia capensis omega virus $T=4$ capsid), `op_id = "2"` and `op_id = "3"` are 5-fold (`C5`, $72^\circ$ / $144^\circ$) pentamer rotations.
- In `6R24.cif` (Saccharomyces cerevisiae Ty3/gypsy retrotransposon capsid), the 5-fold (`C5`, $72^\circ$) pentamer generator is `op_id = "37"`, and `op_id = "4"` is a 3-fold (`C3`, $120^\circ$) rotation.
Therefore, any cross-structure comparison MUST match operators by their invariant rotation order $n_k$, rotation angle $\theta_k$, and axis alignment rather than by raw `operator_id` string.

---

### 4.5 Work Package `L2-E`: `NC` Peptide (`E`) & `Zn` Ion (`F`) Spatial Context and Coordination

#### 4.5.1 Tetrahedral `CCHC` Zinc-Knuckle Coordination Sphere (`C227-X2-C230-X4-H235-X4-C240`)

Within each of the 60 symmetry copies of the C-terminal nucleocapsid (`NC`) fragment (`label_asym_id = E`, residues `224..250`, sequence `TTRCTYCSFRGHTFDNCRKRQKDRQEE`), a single `Zn` ion (`label_asym_id = F`, `auth_asym_id = E`, `auth_seq_id = 301`, ASU coordinates `[373.930, 260.790, 453.921] A`) is tetrahedrally coordinated by three cysteine thiolates and one histidine imidazole nitrogen:

| `_struct_conn.id` | Ligand Residue (`label_asym_id = E`) | `label_seq_id` | `auth_seq_id` | Coordinating Atom | Computed Distance ($\text{\AA}$) | Deposited `_struct_conn` Distance ($\text{\AA}$) |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `metalc1` | `CYS` | 227 | 227 | `SG` | **2.3326** | 2.333 |
| `metalc2` | `CYS` | 230 | 230 | `SG` | **2.3377** | 2.338 |
| `metalc3` | `HIS` | 235 | 235 | `NE2` | **2.1155** | 2.115 |
| `metalc4` | `CYS` | 240 | 240 | `SG` | **2.3482** | 2.348 |

**Tetrahedral Bond Angles at `Zn` (`L1 - ZN - L2`)**:
- `CYS227:SG - ZN - CYS230:SG`: computed `102.77 deg` (deposited `_pdbx_struct_conn_angle`: `102.8 deg`)
- `CYS227:SG - ZN - HIS235:NE2`: computed `105.82 deg` (deposited: `105.8 deg`)
- `CYS230:SG - ZN - HIS235:NE2`: computed `114.40 deg` (deposited: `114.4 deg`)
- `CYS227:SG - ZN - CYS240:SG`: computed `114.03 deg` (deposited: `114.0 deg`)
- `CYS230:SG - ZN - CYS240:SG`: computed `109.39 deg` (deposited: `109.4 deg`)
- `HIS235:NE2 - ZN - CYS240:SG`: computed `110.32 deg` (deposited: `110.3 deg`)
- **Mean tetrahedral angle**: `109.45 deg` (ideal tetrahedral $= 109.47^\circ$); four-coordinate geometry index $\tau_4 = \frac{360^\circ - (\alpha + \beta)}{141^\circ} = \frac{360^\circ - (114.40^\circ + 114.03^\circ)}{141^\circ} = 0.9331$.

#### 4.5.2 `NC-CA` Proximity (`E` to `A, B, C, D`) & `NC-NC` Dimer Proximity (`E[op1] - E[op2]`)

`NC` (`E`) contacts are kept strictly separate from `CA-CA` capsid lattice contacts:
- **`NC-CA` Contacts by Cutoff**:
  - At `4.0 A`: `E-A = 0`, `E-B = 1` (`E[op1]:ARG247 - B[op2]:HIS50`, `3.9941 A`), `E-C = 3` (`E[op1] - C[op1]`, min `3.2124 A` at `GLY234 - THR49`), `E-D = 1` (`E[op1]:ARG247 - D[op1]:TYR91`, `2.6448 A`).
  - At `5.0 A`: `E-A = 0`, `E-B = 1` (`op2`), `E-C = 9` (`op1`), `E-D = 6` (`op1`) (Total = `16` unique `NC-CA` residue pairs per `NC` chain).
  - At `6.0 A`: `E-A = 0`, `E-B = 4` (`1` in `op1`, `3` in `op2`), `E-C = 16` (`op1`), `E-D = 16` (`op1`) (Total = `36` unique `NC-CA` residue pairs per `NC` chain).
- **`NC-NC` Symmetry Dimer Across `op_id = "2"` (`C2` 2-Fold Axis)**:
  - `E[op1]` and `E[op2]` meet beneath the interior of the 2-fold `C-D` / `B-D` hexamer seam (`centroid distance = 15.4794 A`, `Zn[op1]-Zn[op2] distance = 14.8186 A`), forming `21` directed residue-pair contacts (`12` unique unordered residue pairs) at `5.0 A` with minimum heavy-atom distance `3.0194 A` (`E:TYR229 - E:TYR229`).

---

## 5. Mandatory Epistemic Qualifications

1. **`5.0 A` Cutoff is a Geometric Definition, Not a Physical Phase Boundary or Binding Energy**:
   The `5.0 A` heavy-atom distance threshold defines a geometric proximity graph on a static `3.50 A` cryo-EM density model. It is not a thermodynamic binding threshold, covalent bond criterion, or mechanical phase transition. Sensitivity across `3.5 A` to `6.0 A` (`L2-A`) is mandatory for any downstream structural inference.
2. **Both Chain B and Chain D Lack Homotypic Self-Contacts (`B-B = 0` and `D-D = 0`)**:
   B and D both lack homotypic interfaces in the 6TAP residue-contact census (B-B = 0 and D-D = 0 across 3.5..6.0 A; only A-A = 41 and C-C = 7 are homotypic at 5.0 A) alongside `A-C = 0`. In the Level-2 interface graph, `B` separates from `A/C/D` in the exploratory `k=2` partition not by unique absence of homotypic contacts (since `D-D = 0` as well), but by its heterotypic partner-composition and interface-area profile: `B` forms two `B-C` interfaces (`1,272.29` and `1,365.72 A^2`) and two `B-D` interfaces (`1,342.37` and `344.66 A^2`) with only a marginal `48.65 A^2` `A-B` contact, whereas `D` anchors the pentamer-to-hexamer seam via `A-D_op1` (`1,287.13 A^2`, `27` pairs) and `A-D_op36` (`179.98 A^2`, `5` pairs) alongside `C-D_op1` (`1,487.53 A^2`, `34` pairs).
3. **Uniform Subunit Partner Degree (`5.00`) Does Not Imply Uniform Interface Area, Chemistry, or Energetics**:
   Although every `CA` subunit (`A, B, C, D`) contacts exactly $5$ partner subunits at `5.0 A`, the resulting interfaces span a $38.7\times$ range in buried surface area: from `48.65 A^2` (`24.33 A^2` per side, `1` residue pair, `0` H-bond candidates, `0` salt-bridge candidates) at `A-B_C2_peripheral_op2` up to `1,882.00 A^2` (`941.00 A^2` per side, `41` residue pairs, `4` H-bond candidates, `3` salt-bridge candidates) at `A-A_C5_pentamer_op36`. At cutoffs $\le 4.5\text{ \AA}$, `A-B` vanishes (`0` contacts), reducing the partner degree of `A` and `B` from $5$ to $4$.
4. **Zero Transferred Interaction Energies from Unrelated Viral Capsids**:
   No pairwise subunit association energies ($\Delta G$ in kcal/mol or $k_B T$) or assembly rate constants may be copied or rescaled from HIV-1 CA (`3H47`), HBV, Nudaurelia capensis omega virus (`1OHF`), or Ty3/gypsy (`6R24`) into `dArc1` (`6TAP`). Energetic parameterization requires direct experimental measurement on *Drosophila* `dArc1`.

---

## 6. M7 Falsification Design Rules & Verification Disposition

The Level-2 and Level-3 structural interface pipelines enforce the **4 M7 Falsification Design Rules** from the developmental synthesis:

1. **Equifinality requires upstream instrumentation**:
   Macroscopic or downstream summary observables (such as uniform subunit partner degree `5.00`, aggregate capsid diameter, or a `k=2` `60/180` cluster split) are equifinal with respect to distinct microscopic mechanisms. Resolving mechanism requires upstream residue-level, operator-resolved, and multi-metric (`d <= 5.0 A`, Lee-Richards $\Delta\text{SASA}$, Voronoi contact area) instrumentation before any functional claim is made.
2. **Coarse-graining thresholds (`5.0 A`, `22%` foreground) belong to the observer**:
   Discretization cutoffs (`d_cut = 5.0 A` heavy-atom distance, `delta_SASA > 1.0 A^2` burial, or `22%` flow-cytometry foreground gating) are observer-chosen coarse-graining hyper-parameters rather than intrinsic physical phase boundaries. Every structural claim must report sensitivity across the observer cutoff sweep (`3.5 A` to `6.0 A`).
3. **A falsifiable model must beat a competing null model**:
   Structural classifications and mutational hotspot nominations must be benchmarked against explicit competing null models (e.g., sequence-identity / symmetry-equivariant null models, random surface-patch permutation nulls, or area-proportional contact density nulls) with pre-registered experimental falsification readouts.
4. **Validate each link independently before coupling systems**:
   Each stage in the acyclic bridge (`mmCIF coordinate & operator parsing -> independent block-matrix & Voronoi/SASA cross-check -> Level-2 interface atlas -> Level-3 symmetry-normalized interface comparison`) is verified independently with fail-closed receipts before cross-system coupling. Remote publication and CI provenance is recorded explicitly with `"verification_disposition": "remotely_published_on_github_pr6_local_tests_verified_remote_ci_not_run"` and `"publication_and_ci_provenance"`:
   - `"github_pr_url"`: `"https://github.com/aurascoper/flow-cytometry-model-bridge/pull/6"`
   - `"historical_v1_baseline_commit"`: `"955e32dc673d3a874f874ffb8965497db03ca112"`
   - `"level1_v2_repair_commit"`: `"6a09714d03fccf29905c99bddc416785781063d1"`
   - `"level2_atlas_and_crosscheck_commit"`: `"fc5e1a407e375014653b2b72f73b89559c14e2ee"`
   - `"level3_initial_publication_commit"`: `"1abeaad2e44fe25e3e4cd17b7920a51f990fc55e"`
   - `"local_test_execution_summary"`: `"28/28 Python pytest assertions and 2,996/2,996 Julia assertions executed and verified locally"`
   - `"remote_github_actions_ci_status"`: `"No GitHub Actions workflow runs executed on remote; test counts are recorded local verification results, not remote CI reruns"`

