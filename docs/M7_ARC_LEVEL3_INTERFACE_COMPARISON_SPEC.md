# M7 ARC Level-3 Symmetry-Normalized Interface Comparison & 4-Branch Provenance Separation Specification

- **Document ID**: `M7-ARC-L3-INTERFACE-COMPARISON-SPEC-2026-11-07`
- **Schema ID**: `m7-arc-level3-interface-comparison-v1` (`1.3.0`)
- **Verification Disposition**: `"remotely_published_on_github_pr6_local_tests_verified_remote_ci_not_run"`
- **Publication & CI Provenance (`publication_and_ci_provenance`)**:
  - `github_pr_url`: `https://github.com/aurascoper/flow-cytometry-model-bridge/pull/6`
  - `historical_v1_baseline_commit`: `955e32dc673d3a874f874ffb8965497db03ca112`
  - `level1_v2_repair_commit`: `6a09714d03fccf29905c99bddc416785781063d1`
  - `level2_atlas_and_crosscheck_commit`: `fc5e1a407e375014653b2b72f73b89559c14e2ee`
  - `level3_initial_publication_commit`: `1abeaad2e44fe25e3e4cd17b7920a51f990fc55e`
  - `local_test_execution_summary`: `"28/28 Python pytest assertions and 3,059/3,059 Julia assertions (restored from 2,996/2,996) executed and verified locally"`
  - `remote_github_actions_ci_status`: `"No GitHub Actions workflow runs executed on remote; test counts are recorded local verification results, not remote CI reruns"`
- **Primary Structure**: *Drosophila melanogaster* `dArc1` CA Capsid (`PDB 6TAP`, `EMDB EMD-10423`, icosahedral point group `I`, `_em_single_particle_entity.point_symmetry = "I"`, biological assembly `1`, RCSB global symmetry `I`, `60` rotation operators, `240` CA subunits + `60` interior `NC/Zn` segments)
- **Symmetry Provenance Note**: PDB `6TAP` has icosahedral point-group symmetry `I` (`_em_single_particle_entity.point_symmetry = I`, `60` rotation operators). Cubic space-group notation `I432` (No. 211, octahedral) is not applicable to `6TAP`.
- **Upstream Validated Receipts**:
  - `artifacts/m7_arc_level2_validation_receipt.json` (`sha256 = 975fd6a3ada9b3834d82766c0246d6ca8c907f18e65e7bb3faeebcb7bb59432c`)
  - `artifacts/m7_arc_level2_independent_crosscheck_receipt.json` (`sha256 = 0346e16bc5abdf3c657b9781bae4ee0ffa413deea74b4b883b3bb014cfbe9132`)
- **Level-3 Receipt**: `artifacts/m7_arc_level3_interface_comparison_receipt.json` (`sha256 = 17e33810cb298391365c8c589577ae2a45764263e25ac4e5e6434088b67048ee`)

---

## 1. Scope, Provenance Reconciliation & Level-1 Contradiction Repair Summary

Level 3 builds a symmetry-normalized comparative interface atlas directly on top of the 10 symmetry-inequivalent sub-interfaces established in Level 2 (`src/m7_arc_level2_interface_atlas.py`) and cross-checked in Step 3 (`src/m7_arc_level2_independent_crosscheck.py`).

### 1.1 Repaired Level-1 Homotypic-Contact Statement (`:Superseded`)

- **Corrected Census Statement**: `"B and D both lack homotypic interfaces in the 6TAP residue-contact census (B-B = 0 and D-D = 0 across 3.5..6.0 A; only A-A = 41 and C-C = 7 are homotypic at 5.0 A)"`.
- **Disposition of Earlier Level-1 Rationale**: The earlier Level-1 `k=2` narrative claiming that chain `B` was uniquely self-contact-free while `A/C/D` all possessed homotypic contacts is marked `:Superseded`.
- **Level-2/Level-3 Re-evaluation of the `k=2` `B`-vs-`A/C/D` Partition**: Chain `B` separates from `A/C/D` not by unique absence of homotypic contacts (since `D-D = 0` as well), but by its heterotypic partner-composition and interface-area profile:
  - `B` forms two `B-C` interfaces (`B-C_C1_intra_asu_op1` at `1,272.29 A^2`, `26` pairs; and `B-C_C3_inter_asu_hexamer_op33` at `1,365.72 A^2`, `28` pairs) and two `B-D` interfaces (`B-D_C2_dimer_op2` at `1,342.37 A^2`, `28` pairs; and `B-D_C3_inter_asu_vertex_op33` at `344.66 A^2`, `7` pairs) with only a marginal `48.65 A^2` `A-B_C2_peripheral_op2` contact (`1` pair at `4.6059 A`, vanishing at `<= 4.5 A`).
  - Conversely, `D` anchors the pentamer-to-hexamer seam via `A-D_C1_intra_asu_op1` (`1,287.13 A^2`, `27` pairs) and `A-D_C5_inter_asu_seam_op36` (`179.98 A^2`, `5` pairs) alongside `C-D_C1_intra_asu_op1` (`1,487.53 A^2`, `34` pairs).

### 1.2 Computational Implementation Independence Scope & Voronoi-to-BSA Calibration (Step 3)

- **Independence Scope (`COMPUTATIONAL_IMPLEMENTATION_INDEPENDENT_SAME_DEPOSITED_6TAP_COORDINATES`)**: Both `src/m7_arc_level2_interface_atlas.py` and `src/m7_arc_level2_independent_crosscheck.py` operate on the same deposited `PDB 6TAP` atomic coordinates and the same 60 icosahedral assembly operators (`_pdbx_struct_oper_list`). Their independence is strictly **computational and algorithmic** (dense `np.einsum` + Gram-matrix distances + analytical Lee-Richards z-slices + 3D Voronoi radical-plane clipping vs `scipy.spatial.cKDTree` + 100-point Fibonacci Shrake-Rupley SASA), not an independent experimental cryo-EM map reconstruction.
- **Exact Residue-Pair Identity (`Jaccard = 1.0000`)**: All `204` unique asymmetric-unit residue pairs at `5.0 A` (and all census cutoffs `3.5 A` [`50`], `4.0 A` [`79`], `4.5 A` [`136`], `5.0 A` [`204`], `5.5 A` [`284`], and `6.0 A` [`371`]) match 100% on exact canonical residue-pair identities (`symmetric_difference_count = 0`, `jaccard_similarity = 1.0`, `atom_pair_hit_count_mismatches_across_204_pairs = 0`, `max_minimum_distance_abs_diff_A = 4.99e-05 A`).
- **Two-Sided Pairwise BSA vs One-Side Interface Area vs Single-Plane Voronoi Contact Area**:
  - Shrake-Rupley 100-point Fibonacci sphere two-sided pairwise BSA sum across the 10 representative sub-interfaces: `9,570.22 A^2` (`whole_capsid_total_buried_surface_area_A2 = 574,213.20 A^2` across all 60 asymmetric units).
  - Shrake-Rupley one-side interface area sum across the 10 representative sub-interfaces: `4,785.10 A^2` (`whole_capsid_one_side_interface_area_A2 = 287,106.60 A^2` across all 60 asymmetric units).
  - Independent Lee-Richards 40-slice analytical circle-arc contour integration two-sided pairwise BSA sum: `9,400.64 A^2` (`Pearson r = 0.999512` with Shrake-Rupley).
  - Independent Voronoi radical-plane (power-diagram) contact area sum: `4,236.73 A^2` (`Pearson r = 0.999183` with Shrake-Rupley and `r = 0.999763` with Lee-Richards). Because Voronoi contact area integrates a **single dividing radical plane** between contacting subunits whereas two-sided BSA sums solvent-accessible surface lost on **both** subunits (`BSA_two_sided = 2 * BSA_one_side`), the OLS linear calibration across the 10 sub-interfaces is:
    - Against two-sided Shrake-Rupley BSA: `voronoi_A2 = 0.444876 * sr_pairwise_bsa_A2 - 2.0833` (`R^2 = 0.998367`, `RMSE = 11.1209 A^2`, `max |residual| = 24.1430 A^2`).
    - Against one-side Shrake-Rupley BSA: `voronoi_A2 = 0.889755 * sr_one_side_bsa_A2 - 2.0838` (`R^2 = 0.998368`, `RMSE = 11.1203 A^2`, `max |residual| = 24.1445 A^2`).

### 1.3 Exact mmCIF `_entity`, `_entity_poly`, `_struct_asym`, and `_atom_site.label_entity_id` Reconciliation for Chains `A..F`

Direct inspection of `structures/6TAP.cif` establishes the exact mmCIF provenance of all six asymmetric-unit strands (`A..F`):
- **`_entity.id = "1"`**: `polymer`, `polypeptide(L)`, `'Activity-regulated cytoskeleton associated protein 1'` (`dArc1`, UniProt `Q7K1U0`, `254` aa, `_entity.pdbx_number_of_molecules = 5`, `_entity_poly.pdbx_strand_id = A,B,C,D,E`).
- **`_entity.id = "2"`**: `non-polymer`, `'ZINC ION'` (`comp_id = ZN`, `_entity.pdbx_number_of_molecules = 1`, `_pdbx_nonpoly_scheme.asym_id = F`, `_pdbx_nonpoly_scheme.pdb_strand_id = E`, `pdb_seq_num = 301`, `auth_seq_num = 251`).
- Chain `E` (`NC` zinc knuckle, residues `224..250`) is a resolved C-terminal segment of the **same** deposited `dArc1` polypeptide Entity `1` (NOT an independently deposited polymer `"Entity 2"`), and the `ZN` ion (`label_asym_id = F`, `auth_asym_id = E`) is non-polymer ligand Entity `2` (NOT `"Entity 3"`). Functional roles are distinguished via `functional_component_id` and `biological_domain_or_fragment_designation`:

| `label_asym_id` | `auth_asym_id` | `label_entity_id` | `struct_asym_entity_id` | `functional_component_id` | `biological_domain_or_fragment_designation` | Resolved `_atom_site` Range & Atom Count |
| :---: | :---: | :---: | :---: | :--- | :--- | :--- |
| `A` | `A` | `"1"` | `"1"` | `"CA_CAPSID_SHELL"` | `"CA_capsid_shell_domain_pentamer_A_41_205"` | `42..205` (`164` aa, `1,343` atoms) |
| `B` | `B` | `"1"` | `"1"` | `"CA_CAPSID_SHELL"` | `"CA_capsid_shell_domain_hexamer_B_43_205"` | `41..204` (`164` aa, `1,339` atoms) |
| `C` | `C` | `"1"` | `"1"` | `"CA_CAPSID_SHELL"` | `"CA_capsid_shell_domain_hexamer_C_43_206"` | `41..205` (`165` aa, `1,347` atoms) |
| `D` | `D` | `"1"` | `"1"` | `"CA_CAPSID_SHELL"` | `"CA_capsid_shell_domain_hexamer_D_41_206"` | `41..206` (`166` aa, `1,356` atoms) |
| `E` | `E` | `"1"` | `"1"` | `"NC_ZINC_KNUCKLE_FRAGMENT"` | `"NC_zinc_knuckle_C_terminal_fragment_224_250_of_entity_1"` | `224..250` (`27` aa, `233` atoms) |
| `F` | `E` | `"2"` | `"2"` | `"ZN_COORDINATION_LIGAND"` | `"ZN_tetrahedral_CCHC_ligand_entity_2"` | `_atom_site.auth_seq_id 301` / `pdb_seq_num 301` / `auth_seq_num 251` (`1` `ZN` atom) |

---

## 2. Symmetry-Normalized Residue Fingerprints & Interface-Fingerprint Families Across the 10 Sub-Interfaces

Because chains `A` (`42..205`, `164` aa), `B` (`41..204`, `164` aa), `C` (`41..205`, `165` aa), and `D` (`41..206`, `166` aa; `659` total CA residues across the ASU) are quasi-equivalent copies of the same `dArc1` CA sequence (`label_entity_id = "1"`), Level 3 projects every interfacial residue onto its canonical sequence coordinate (`label_seq_id == auth_seq_id`), stratified into:
- **NTD (`seq_id 41..125`)**: NTD stalk (`41..54`) + NTD four-helix bundle (`55..125`)
- **CTD (`seq_id 126..206`)**: Interdomain hinge (`126..132`) + CTD five-helix lobe/tail (`133..206`)

> **Terminology Note (`interface_fingerprint_family`)**: Throughout Level 3, Families I, II, and III (`Family_III_Residual_Distal_CTD_Hierarchy`) are strictly referred to as **interface-fingerprint families** / **interface-fingerprint groupings** (`interface_fingerprint_family`), i.e., empirical groupings of sub-interfaces defined by homologous sequence-position similarity and NTD/CTD contact architecture in `6TAP`, **not** established biological or functional modules, and are explicitly stratified by their differential threshold robustness (`family_robustness_stratification`).

### 2.1 Global Cross-Chain Homologous Residue Conservation

- **Total unique homologous sequence positions (`seq_id`) participating in at least one `<= 5.0 A` contact across `A, B, C, D`**: `83` unique sequence positions (out of `166` resolved sequence positions `41..206`).
- **Per-chain interfacial residue counts**:
  - Chain `A`: `62` interfacial residues
  - Chain `B`: `53` interfacial residues
  - Chain `C`: `60` interfacial residues
  - Chain `D`: `60` interfacial residues
- **Conserved interfacial sequence positions active across all 4 chains (`A, B, C, D`)** (`35` positions):
  `46, 56, 70, 74, 88, 93, 95, 96, 97, 100, 101, 117, 118, 125, 126, 129, 130, 133, 144, 145, 148, 151, 152, 161, 162, 166, 169, 170, 172, 173, 174, 185, 195, 199, 204`
- **Conserved interfacial sequence positions active across `>= 3` chains** (`51` positions):
  `44, 46, 56, 59, 63, 67, 70, 71, 73, 74, 75, 84, 87, 88, 93, 95, 96, 97, 100, 101, 104, 117, 118, 121, 125, 126, 129, 130, 133, 144, 145, 148, 151, 152, 161, 162, 166, 169, 170, 172, 173, 174, 180, 181, 185, 192, 195, 199, 200, 203, 204`

### 2.2 Sub-Interface Residue & Interface-Fingerprint Family Table

| Sub-Interface Class | Interface-Fingerprint Family (`interface_fingerprint_family`) | Subfamily Tier (`interface_fingerprint_subfamily_tier`) | Chain Pair | Sym (`op`) | Pairs (`5.0 A`) | Chain Res (`NTD` / `CTD`) | Homologous `seq_id` Count | NTD `seq_id` Set (`41..125`) | CTD `seq_id` Set (`126..206`) |
| :--- | :--- | :--- | :---: | :---: | :---: | :---: | :---: | :--- | :--- |
| `A-A_C5_pentamer_op36` | `Family_I_NTD_CTD_Capsomere_Seam` | `Family_I_NTD_CTD_Architectural_Seam` | `A-A` | `C5` (`1,36`) | `41` | `44` (`30` / `14`) | `44` | `44, 46, 50, 56, 59, 63, 67, 70, 71, 73, 74, 75, 84, 85, 87, 88, 90, 91, 92, 93, 95, 96, 97, 100, 101, 104, 117, 118, 121, 124` | `126, 142, 144, 145, 148, 151, 152, 192, 195, 196, 199, 202, 203, 204` |
| `C-D_C1_intra_asu_op1` | `Family_I_NTD_CTD_Capsomere_Seam` | `Family_I_NTD_CTD_Architectural_Seam` | `C-D` | `C1` (`1,1`) | `34` | `42` (`28` / `14`) | `42` | `41, 42, 44, 46, 49, 50, 56, 63, 70, 71, 73, 74, 75, 81, 84, 87, 88, 93, 95, 96, 97, 99, 100, 101, 104, 117, 118, 121` | `127, 144, 145, 148, 149, 151, 152, 155, 156, 164, 192, 195, 199, 203` |
| `B-C_C3_inter_asu_hexamer_op33` | `Family_II_Pure_CTD_Quasi_Equivalent_Dimer` | `Family_II_Coherent_CTD_Dimer_Core` | `B-C` | `C3` (`1,33`) | `28` | `31` (`2` / `29`) | `17` | `125` | `126, 129, 130, 133, 161, 162, 165, 166, 167, 169, 170, 172, 173, 174, 180, 185` |
| `B-D_C2_dimer_op2` | `Family_I_NTD_CTD_Capsomere_Seam` | `Family_I_NTD_CTD_Architectural_Seam` | `B-D` | `C2` (`1,2`) | `28` | `33` (`25` / `8`) | `33` | `44, 46, 56, 59, 63, 67, 70, 71, 73, 74, 75, 84, 87, 88, 90, 91, 92, 93, 95, 96, 97, 100, 101, 117, 118` | `142, 144, 145, 148, 151, 152, 195, 199` |
| `A-D_C1_intra_asu_op1` | `Family_II_Pure_CTD_Quasi_Equivalent_Dimer` | `Family_II_Coherent_CTD_Dimer_Core` | `A-D` | `C1` (`1,1`) | `27` | `30` (`2` / `28`) | `16` | `125` | `126, 129, 130, 133, 136, 161, 162, 166, 169, 170, 172, 173, 174, 180, 185` |
| `B-C_C1_intra_asu_op1` | `Family_I_NTD_CTD_Capsomere_Seam` | `Family_I_NTD_CTD_Architectural_Seam` | `B-C` | `C1` (`1,1`) | `26` | `34` (`21` / `13`) | `34` | `42, 46, 56, 58, 59, 66, 67, 70, 74, 88, 93, 95, 96, 97, 99, 100, 101, 104, 117, 118, 121` | `126, 127, 144, 145, 148, 151, 152, 155, 164, 192, 195, 199, 203` |
| `C-C_C3_trimer_op33` | `Family_III_Residual_Distal_CTD_Hierarchy` | `Family_III_Tier1_Coherent_C3_Vertex_Core` | `C-C` | `C3` (`1,33`) | `7` | `9` (`0` / `9`) | `9` | *(none)* | `177, 180, 181, 182, 186, 197, 200, 204, 205` |
| `B-D_C3_inter_asu_vertex_op33` | `Family_III_Residual_Distal_CTD_Hierarchy` | `Family_III_Tier1_Coherent_C3_Vertex_Core` | `B-D` | `C3` (`1,33`) | `7` | `9` (`0` / `9`) | `9` | *(none)* | `177, 180, 181, 182, 184, 200, 201, 204, 206` |
| `A-D_C5_inter_asu_seam_op36` | `Family_III_Residual_Distal_CTD_Hierarchy` | `Family_III_Tier2_C5_Inter_ASU_Seam` | `A-D` | `C5` (`1,36`) | `5` | `6` (`0` / `6`) | `6` | *(none)* | `184, 185, 186, 200, 204, 205` |
| `A-B_C2_peripheral_op2` | `Family_III_Residual_Distal_CTD_Hierarchy` | `Family_III_Tier3_Peripheral_C2_Singleton` | `A-B` | `C2` (`1,2`) | `1` | `2` (`0` / `2`) | `2` | *(none)* | `181, 204` |

### 2.3 $10 \times 10$ Pairwise Jaccard Similarity Matrix (Homologous `seq_id` Positions)

| Sub-Interface | `A-A_op36` | `C-D_op1` | `B-C_op33` | `B-D_op2` | `A-D_op1` | `B-C_op1` | `C-C_op33` | `B-D_op33` | `A-D_op36` | `A-B_op2` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `A-A_C5_pentamer_op36` | `1.0000` | `0.5926` | `0.0167` | `0.7500` | `0.0169` | `0.5294` | `0.0192` | `0.0192` | `0.0204` | `0.0222` |
| `C-D_C1_intra_asu_op1` | `0.5926` | `1.0000` | `0.0000` | `0.5625` | `0.0000` | `0.6170` | `0.0000` | `0.0000` | `0.0000` | `0.0000` |
| `B-C_C3_inter_asu_hexamer_op33` | `0.0167` | `0.0000` | `1.0000` | `0.0000` | `0.8333` | `0.0200` | `0.0400` | `0.0400` | `0.0455` | `0.0000` |
| `B-D_C2_dimer_op2` | `0.7500` | `0.5625` | `0.0000` | `1.0000` | `0.0000` | `0.4889` | `0.0000` | `0.0000` | `0.0000` | `0.0000` |
| `A-D_C1_intra_asu_op1` | `0.0169` | `0.0000` | `0.8333` | `0.0000` | `1.0000` | `0.0204` | `0.0417` | `0.0417` | `0.0476` | `0.0000` |
| `B-C_C1_intra_asu_op1` | `0.5294` | `0.6170` | `0.0200` | `0.4889` | `0.0204` | `1.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` |
| `C-C_C3_trimer_op33` | `0.0192` | `0.0000` | `0.0400` | `0.0000` | `0.0417` | `0.0000` | `1.0000` | `0.5000` | `0.3636` | `0.2222` |
| `B-D_C3_inter_asu_vertex_op33` | `0.0192` | `0.0000` | `0.0400` | `0.0000` | `0.0417` | `0.0000` | `0.5000` | `1.0000` | `0.2500` | `0.2222` |
| `A-D_C5_inter_asu_seam_op36` | `0.0204` | `0.0000` | `0.0455` | `0.0000` | `0.0476` | `0.0000` | `0.3636` | `0.2500` | `1.0000` | `0.1429` |
| `A-B_C2_peripheral_op2` | `0.0222` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.2222` | `0.2222` | `0.1429` | `1.0000` |

### 2.4 $10 \times 10$ Pairwise Jaccard Similarity Matrix (Canonical Homologous Pairs `[min_seq, max_seq]`)

| Sub-Interface | `A-A_op36` | `C-D_op1` | `B-C_op33` | `B-D_op2` | `A-D_op1` | `B-C_op1` | `C-C_op33` | `B-D_op33` | `A-D_op36` | `A-B_op2` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `A-A_C5_pentamer_op36` | `1.0000` | `0.2097` | `0.0000` | `0.6429` | `0.0000` | `0.2407` | `0.0000` | `0.0000` | `0.0000` | `0.0000` |
| `C-D_C1_intra_asu_op1` | `0.2097` | `1.0000` | `0.0000` | `0.1698` | `0.0000` | `0.4286` | `0.0000` | `0.0000` | `0.0000` | `0.0000` |
| `B-C_C3_inter_asu_hexamer_op33` | `0.0000` | `0.0000` | `1.0000` | `0.0000` | `0.6364` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` |
| `B-D_C2_dimer_op2` | `0.6429` | `0.1698` | `0.0000` | `1.0000` | `0.0000` | `0.2000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` |
| `A-D_C1_intra_asu_op1` | `0.0000` | `0.0000` | `0.6364` | `0.0000` | `1.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` |
| `B-C_C1_intra_asu_op1` | `0.2407` | `0.4286` | `0.0000` | `0.2000` | `0.0000` | `1.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` |
| `C-C_C3_trimer_op33` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `1.0000` | `0.4000` | `0.0909` | `0.1429` |
| `B-D_C3_inter_asu_vertex_op33` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.4000` | `1.0000` | `0.0000` | `0.1429` |
| `A-D_C5_inter_asu_seam_op36` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0909` | `0.0000` | `1.0000` | `0.0000` |
| `A-B_C2_peripheral_op2` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.0000` | `0.1429` | `0.1429` | `0.0000` | `1.0000` |

### 2.5 Differential Robustness of the Three Interface-Fingerprint Groupings (`family_robustness_stratification`)

The three empirical interface-fingerprint groupings in `PDB 6TAP` exhibit **sharply unequal threshold robustness** across homologous sequence-position (`seq_id`) vs canonical residue-pair comparisons and single-linkage vs complete-linkage clustering:

| Grouping (`interface_fingerprint_family`) | Robustness Rank & Tier | Members | Homologous `seq_id` Jaccard Grid Interval (Single / Complete) | Canonical Residue-Pair Jaccard Grid Interval (Single / Complete) | Structural Qualification & Hierarchical Resolution |
| :--- | :--- | :--- | :---: | :---: | :--- |
| **`Family II` (`Family_II_Pure_CTD_Quasi_Equivalent_Dimer`)** | **Rank 1**: `STRONGEST_SHARED_CTD_INTERFACE_FINGERPRINT_FAMILY` | `A-D_C1_intra_asu_op1`, `B-C_C3_inter_asu_hexamer_op33` | `[0.15, 0.70]` / `[0.15, 0.70]` (continuous `(0.0476, 0.8333]`) | `[0.15, 0.60]` / `[0.15, 0.60]` (continuous `(0.0000, 0.6364]`) | **Strongest shared CTD interface-fingerprint family**: survives broad thresholds across both homologous sequence-position (`15/18` shared `seq_ids`, `Jaccard = 0.8333`) and canonical residue-pair (`14` shared canonical pairs, `Jaccard = 0.6364`) comparisons with `0.0000` canonical-pair overlap against all other sub-interfaces. |
| **`Family I` (`Family_I_NTD_CTD_Capsomere_Seam`)** | **Rank 2**: `BROADER_NTD_CTD_ARCHITECTURAL_SEAM_MOTIF_THRESHOLD_SENSITIVE_AT_EXACT_PAIR_LEVEL` | `A-A_C5_pentamer_op36`, `C-D_C1_intra_asu_op1`, `B-C_C1_intra_asu_op1`, `B-D_C2_dimer_op2` | `[0.15, 0.50]` / `[0.15, 0.45]` (continuous `(0.0222, 0.5926]` / `(0.0222, 0.4889]`) | `[0.15, 0.20]` / `[0.15, 0.15]` (continuous `(0.0000, 0.2407]` / `(0.0000, 0.1698]`) | **Broader NTD/CTD architectural seam motif**: strongly supported at the homologous sequence-position / domain-architecture level (`[0.15, 0.50]` single linkage, `[0.15, 0.45]` complete linkage Jaccard), whereas exact canonical residue-pair correspondence is substantially more threshold-sensitive (`[0.15, 0.20]` single linkage and `0.15` complete linkage Jaccard) because quasi-equivalent rotational registers shift individual partner pairings across the shared NTD/CTD seam. |
| **`Family III` (`Family_III_Residual_Distal_CTD_Hierarchy`)** | **Rank 3**: `RESIDUAL_DISTAL_CTD_INTERFACE_HIERARCHY_NOT_COHESIVE_FOUR_MEMBER_FAMILY` | `C-C_C3_trimer_op33`, `B-D_C3_inter_asu_vertex_op33`, `A-D_C5_inter_asu_seam_op36`, `A-B_C2_peripheral_op2` | `[0.15, 0.20]` single / **never** (`[]`) at `tau >= 0.15` complete | **Never** (`[]`) at `tau >= 0.15` single or complete | **Residual distal CTD interface hierarchy (NOT one cohesive four-member family)**: never forms a four-member cluster at any tested complete-linkage threshold $\tau \ge 0.15$, resolving into: (1) a coherent two-member **`C3` vertex core** (`{C-C_op33, B-D_op33}`, Jaccard `0.5000` on `seq_id`, `0.4000` on canonical pairs), (2) the **`C5` inter-ASU seam** (`A-D_op36`, detaching above $\tau = 0.2500$ complete / $0.3636$ single Jaccard on `seq_id` and $\tau = 0.0000$ complete / $0.0909$ single on canonical pairs), and (3) the **peripheral `C2` singleton** (`A-B_op2`, detaching above $\tau = 0.1429$ complete / $0.2222$ single Jaccard on `seq_id` and $\tau = 0.0000$ complete / $0.1429$ single on canonical pairs). |

- **Same-Chain-Pair Sub-Interfaces Use Completely Distinct Surface Patches**:
  - **`B-C_C1_intra_asu_op1` vs `B-C_C3_inter_asu_hexamer_op33`**: Canonical residue-pair Jaccard is **`0.0000`** (`0` shared pairs) and `seq_id` Jaccard is **`0.0200`** (sharing only hinge residue `B:TYR126`; on chain `C`, overlap is `0.0000`).
  - **`A-D_C1_intra_asu_op1` vs `A-D_C5_inter_asu_seam_op36`**: Canonical residue-pair Jaccard is **`0.0000`** and `seq_id` Jaccard is **`0.0476`** (sharing only `D:ARG185`; on chain `A`, overlap is `0.0000`).
  - **`B-D_C2_dimer_op2` vs `B-D_C3_inter_asu_vertex_op33`**: Both `seq_id` Jaccard and canonical pair Jaccard are **`0.0000`** (`0` shared residues on `B` and `0` shared residues on `D`).

### 2.6 Alternative Residue-Similarity Metrics & Clustering Threshold-Sweep Sensitivity (`alternative_similarity_metrics_and_threshold_sensitivity`)

To verify how the three groupings behave across alternative similarity definitions and sub-interface size disparity (`1..41` residue pairs), `src/m7_arc_level3_interface_comparison.py` evaluates $10 \times 10$ similarity matrices across **4 set-similarity metrics** on both homologous `seq_id` sets and canonical homologous residue-pair sets, plus **contact-degree-weighted cosine similarity** across the `83` homologous sequence positions:
1. **`jaccard`**: $|A \cap B| / |A \cup B|$
2. **`sorensen_dice`**: $2|A \cap B| / (|A| + |B|)$
3. **`szymkiewicz_simpson_overlap`**: $|A \cap B| / \min(|A|, |B|)$
4. **`cosine_ochiai`**: $|A \cap B| / \sqrt{|A||B|}$
5. **`contact_degree_weighted_cosine`**: $\sum_s w_u(s) w_v(s) / (\|w_u\|_2 \|w_v\|_2)$ across the `83` homologous sequence positions.

Single-linkage and complete-linkage clustering sweeps across $\tau \in \{0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.60, 0.70\}$ establish the exact threshold intervals $[\tau_{\min}, \tau_{\max}]$ where Family II and Family I remain invariant and how `Family_III_Residual_Distal_CTD_Hierarchy` resolves into its three tiers:

| Representation & Metric | Family I `{A-A_op36, C-D_op1, B-C_op1, B-D_op2}` Single-Linkage Grid $[\tau_{\min}, \tau_{\max}]$ (Continuous) | Family I Complete-Linkage Grid $[\tau_{\min}, \tau_{\max}]$ (Continuous) | Family II `{A-D_op1, B-C_op33}` Single & Complete Grid $[\tau_{\min}, \tau_{\max}]$ (Continuous) | Family III Tier 3 `A-B_op2` Detachment Threshold (Complete / Single) | Family III Tier 2 `A-D_op36` Detachment Threshold (Complete / Single) | Family III Tier 1 `C3` Core `{C-C_op33, B-D_op33}` Similarity |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **Homologous `seq_id`: `jaccard`** | `[0.15, 0.50]` (`(0.0222, 0.5926]`) | `[0.15, 0.45]` (`(0.0222, 0.4889]`) | `[0.15, 0.70]` (`(0.0476, 0.8333]`) | `0.1429` / `0.2222` | `0.2500` / `0.3636` | `0.5000` |
| **Homologous `seq_id`: `sorensen_dice`** | `[0.15, 0.70]` (`(0.0435, 0.7442]`) | `[0.15, 0.60]` (`(0.0435, 0.6567]`) | `[0.15, 0.70]` (`(0.0909, 0.9091]`) | `0.2500` / `0.3636` | `0.4000` / `0.5333` | `0.6667` |
| **Homologous `seq_id`: `szymkiewicz_simpson_overlap`** | `[0.60, 0.70]` (`(0.5000, 0.8182]`) | `[0.15, 0.60]` (`(0.0000, 0.6667]`) | `[0.20, 0.70]` (`(0.1667, 0.9375]`) | `0.5000` / `1.0000` | `0.5000` / `0.6667` | `0.6667` |
| **Homologous `seq_id`: `cosine_ochiai`** | `[0.15, 0.70]` (`(0.1066, 0.7444]`) | `[0.15, 0.60]` (`(0.1066, 0.6568]`) | `[0.15, 0.70]` (`(0.1021, 0.9095]`) | `0.2887` / `0.4714` | `0.4082` / `0.5443` | `0.6667` |
| **Homologous `seq_id`: `contact_degree_weighted_cosine`** | `[0.15, 0.70]` (`(0.0558, 0.7036]`) | `[0.15, 0.60]` (`(0.0558, 0.6221]`) | `[0.15, 0.70]` (`(0.0558, 0.9571]`) | `0.4743` / `0.8750` | `0.4743` / `0.6339` | `0.8352` |
| **Canonical Pair: `jaccard`** | `[0.15, 0.20]` (`(0.0000, 0.2407]`) | `[0.15, 0.15]` (`(0.0000, 0.1698]`) | `[0.15, 0.60]` (`(0.0000, 0.6364]`) | `0.0000` / `0.1429` | `0.0000` / `0.0909` | `0.4000` |
| **Canonical Pair: `sorensen_dice`** | `[0.15, 0.35]` (`(0.0000, 0.3881]`) | `[0.15, 0.25]` (`(0.0000, 0.2903]`) | `[0.15, 0.70]` (`(0.0000, 0.7778]`) | `0.0000` / `0.2500` | `0.0000` / `0.1667` | `0.5714` |
| **Canonical Pair: `szymkiewicz_simpson_overlap`** | `[0.15, 0.50]` (`(0.0000, 0.5000]`) | `[0.15, 0.30]` (`(0.0000, 0.3214]`) | `[0.15, 0.70]` (`(0.0000, 0.8235]`) | `0.0000` / `1.0000` | `0.0000` / `0.2000` | `0.5714` |
| **Canonical Pair: `cosine_ochiai`** | `[0.15, 0.35]` (`(0.0000, 0.3982]`) | `[0.15, 0.25]` (`(0.0000, 0.2917]`) | `[0.15, 0.70]` (`(0.0000, 0.7790]`) | `0.0000` / `0.3780` | `0.0000` / `0.1690` | `0.5714` |

---

## 3. Normalized Interface Chemistry & Area Profiles

### 3.1 Area & Contact Density Comparison Table Across the 10 Sub-Interfaces

| Sub-Interface Class | Sym (`op`) | Shrake-Rupley Two-Sided BSA (`A^2`) | One-Side BSA (`A^2`) | Lee-Richards BSA (`A^2`) | Voronoi Area (`A^2`) | Whole-Capsid Two-Sided BSA (`A^2`) | Whole-Capsid One-Side BSA (`A^2`) | Pairs (`5 A`) | Rep Atom Hits (`5 A`) | Atom Hits / `100 A^2` | Pairs / `1,000 A^2` |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `A-A_C5_pentamer_op36` | `C5` (`1,36`) | `1,882.00` | `941.00` | `1,843.14` | `837.98` | `112,920.00` | `56,460.00` | `41` | `218.0` | `11.58` | `21.79` |
| `C-D_C1_intra_asu_op1` | `C1` (`1,1`) | `1,487.53` | `743.76` | `1,472.78` | `660.10` | `89,251.80` | `44,625.60` | `34` | `143.0` | `9.61` | `22.86` |
| `B-C_C3_inter_asu_hexamer_op33` | `C3` (`1,33`) | `1,365.72` | `682.86` | `1,322.38` | `581.35` | `81,943.20` | `40,971.60` | `28` | `252.0` | `18.45` | `20.50` |
| `B-D_C2_dimer_op2` | `C2` (`1,2`) | `1,342.37` | `671.18` | `1,345.08` | `611.62` | `80,542.20` | `40,270.80` | `28` | `106.0` | `7.90` | `20.86` |
| `A-D_C1_intra_asu_op1` | `C1` (`1,1`) | `1,287.13` | `643.56` | `1,254.46` | `569.95` | `77,227.80` | `38,613.60` | `27` | `185.0` | `14.37` | `20.98` |
| `B-C_C1_intra_asu_op1` | `C1` (`1,1`) | `1,272.29` | `636.15` | `1,278.56` | `572.26` | `76,337.40` | `38,169.00` | `26` | `121.0` | `9.51` | `20.44` |
| `C-C_C3_trimer_op33` | `C3` (`1,33`) | `359.89` | `179.94` | `325.12` | `146.43` | `21,593.40` | `10,796.40` | `7` | `25.5` | `7.09` | `19.45` |
| `B-D_C3_inter_asu_vertex_op33` | `C3` (`1,33`) | `344.66` | `172.33` | `319.29` | `147.40` | `20,679.60` | `10,339.80` | `7` | `27.0` | `7.83` | `20.31` |
| `A-D_C5_inter_asu_seam_op36` | `C5` (`1,36`) | `179.98` | `89.99` | `171.41` | `77.67` | `10,798.80` | `5,399.40` | `5` | `26.0` | `14.45` | `27.78` |
| `A-B_C2_peripheral_op2` | `C2` (`1,2`) | `48.65` | `24.33` | `68.42` | `31.97` | `2,919.00` | `1,459.80` | `1` | `3.0` | `6.17` | `20.55` |
| **Total / Sum** | **10 Classes** | **`9,570.22`** | **`4,785.10`** | **`9,400.64`** | **`4,236.73`** | **`574,213.20`** | **`287,106.60`** | **`204`** | **`1,106.5`** | **`11.56` (mean)** | **`21.32` (mean)** |

### 3.2 Candidate Chemistry Density Profile (Per `1,000 A^2` Shrake-Rupley BSA)

All interaction chemistry counts are explicitly qualified as **geometric candidate counts** (`hydrophobic_candidate_pairs` at `<= 4.5 A` side-chain C/S; `pi_stacking_candidate_pairs` at `<= 7.0 A` ring centroid with plane angle; `hbond_heavy_atom_distance_candidate_pairs` at `<= 3.5 A` N/O heavy-atom distance without explicit H-atom donor-H-acceptor angle; `salt_bridge_heavy_atom_distance_candidate_pairs` at `<= 4.0 A` cationic-N to anionic-O distance):

| Sub-Interface Class | Hydrophobic Cand. (`count` / `per 1k A^2`) | Pi-Stacking Cand. (`count` / `per 1k A^2`) | H-Bond Heavy-Atom Cand. (`count` / `per 1k A^2`) | Salt-Bridge Heavy-Atom Cand. (`count` / `per 1k A^2`) | Domain Mode (`NTD-NTD` / `NTD-CTD` / `CTD-CTD`) |
| :--- | :---: | :---: | :---: | :---: | :---: |
| `A-A_C5_pentamer_op36` | `2` / `1.06` | `0` / `0.00` | `4` / `2.13` | `3` / `1.59` | Composite (`20` / `20` / `1`) |
| `C-D_C1_intra_asu_op1` | `1` / `0.67` | `0` / `0.00` | `3` / `2.02` | `3` / `2.02` | Composite (`12` / `21` / `1`) |
| `B-C_C3_inter_asu_hexamer_op33` | `6` / `4.39` | `4` / `2.93` | `5` / `3.66` | `0` / `0.00` | CTD-Dimer (`0` / `2` / `26`) |
| `B-D_C2_dimer_op2` | `3` / `2.23` | `0` / `0.00` | `0` / `0.00` | `1` / `0.74` | Composite (`17` / `11` / `0`) |
| `A-D_C1_intra_asu_op1` | `7` / `5.44` | `4` / `3.11` | `1` / `0.78` | `0` / `0.00` | CTD-Dimer (`0` / `2` / `25`) |
| `B-C_C1_intra_asu_op1` | `1` / `0.79` | `0` / `0.00` | `3` / `2.36` | `2` / `1.57` | Composite (`7` / `17` / `2`) |
| `C-C_C3_trimer_op33` | `0` / `0.00` | `1` / `2.78` | `1` / `2.78` | `0` / `0.00` | Pure CTD (`0` / `0` / `7`) |
| `B-D_C3_inter_asu_vertex_op33` | `0` / `0.00` | `1` / `2.90` | `0` / `0.00` | `0` / `0.00` | Pure CTD (`0` / `0` / `7`) |
| `A-D_C5_inter_asu_seam_op36` | `0` / `0.00` | `0` / `0.00` | `2` / `11.11` | `0` / `0.00` | Pure CTD (`0` / `0` / `5`) |
| `A-B_C2_peripheral_op2` | `0` / `0.00` | `0` / `0.00` | `0` / `0.00` | `0` / `0.00` | Pure CTD (`0` / `0` / `1`) |
| **ASU Total** | **`20` (`2.09 / 1k A^2`)** | **`10` (`1.04 / 1k A^2`)** | **`19` (`1.99 / 1k A^2`)** | **`9` (`0.94 / 1k A^2`)** | **`56` / `73` / `75` (`204`)** |

### 3.3 Competing Null Model Evaluations, Non-Independence Caveats & Degree/Interface-Preserving Permutation Nulls (M7 Falsification Design Rule 3)

Per M7 Falsification Design Rule 3 (*"A falsifiable model must beat a competing null model"*), Level 3 evaluates both **model-based exploratory screening nulls** and **degree- and interface-membership-preserving permutation nulls** (`"inference_tier": "MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY"`, epistemic status `:Incomplete` until wet-lab mutagenesis).

#### 3.3.1 Non-Independence Caveats (`non_independence_caveats`)

The `204` symmetry-unique residue-pair contacts at `d <= 5.0 A` are **not statistically independent i.i.d. draws**:
1. **Vertex-degree recurrence**: The `204` contact edges connect `235` chain-specific `(chain_id, seq_id)` nodes and `83` homologous sequence positions (`38` NTD, `45` CTD) with heterogeneous contact degrees (up to degree `13`).
2. **Sub-interface rigid-body coupling**: Multiple contacts belong to the same sub-interface (partitioned across `10` sub-interfaces ranging from `1` to `41` residue pairs) and are coupled by rigid-body subunit placement.
3. **Exact icosahedral symmetry copies**: The `60` asymmetric units in `PDB 6TAP` (`point_symmetry = "I"`) are exact geometric copies generated by `_pdbx_struct_oper_list` rather than independent biological replicates.
4. **Exploratory screening scope**: Consequently, both the multinomial domain-segregation chi-square ($\chi^2 = 22.1774$, nominal $p_{\text{iid}} = 1.53 \times 10^{-5}$) and the area-proportional chemistry null (`2.35x` hydrophobic, `2.89x` pi-stacking in Family II) are explicitly labeled `"inference_tier": "MODEL_BASED_EXPLORATORY_SCREENING_NULL_NOT_CONFIRMATORY"`.

#### 3.3.2 Model-Based Exploratory Screening Null Comparisons

1. **Multinomial Domain-Segregation Screening Null (`domain_segregation_vs_random_mixing_null`)**:
   - Across the `659` resolved CA residues in the `6TAP` ASU, `336/659` (`50.99%`) belong to the NTD (`seq_id <= 125`) and `323/659` (`49.01%`) belong to the CTD (`seq_id >= 126`).
   - Under a nominal i.i.d. random domain-pairing baseline across `204` contacts, expected counts are `NTD-NTD = 53.03`, `NTD-CTD = 101.96`, and `CTD-CTD = 49.01`.
   - Observed counts across the `204` contacts are `NTD-NTD = 56`, `NTD-CTD = 73`, and `CTD-CTD = 75` ($\chi^2 = 22.1774$, $\text{df} = 2$, nominal $p_{\text{iid}} = 1.53 \times 10^{-5}$): four distal interfaces in `Family_III_Residual_Distal_CTD_Hierarchy` (`C-C_op33`, `B-D_op33`, `A-D_op36`, `A-B_op2`) are `100%` pure `CTD-CTD` (`20/20` pairs) and the two CTD dimer locks in Family II (`A-D_op1` and `B-C_op33`) are `92.7%` `CTD-CTD` (`51/55` pairs).
2. **Area-Proportional Chemistry Screening Null (`ctd_dimer_hydrophobic_aromatic_vs_ntd_seam_salt_bridge_null`)**:
   - Under a uniform area-proportional baseline, Family II (`A-D_op1` + `B-C_op33`, combined `BSA = 2,652.85 A^2` = `27.72%` of total representative BSA) would carry `5.54` hydrophobic candidates, `2.77` pi-stacking candidates, and `2.49` salt-bridge candidates.
   - Observed Family II counts are **`13/20` (`65.0%`, `2.35x` enrichment) of all hydrophobic candidates**, **`8/10` (`80.0%`, `2.89x` enrichment) of all pi-stacking candidates**, and **`0/9` (`0.0%`) of salt-bridge candidates**.
   - Conversely, **`9/9` (`100.0%`, `1.60x` enrichment) of salt-bridge atom-pair candidates** (`6/6` salt-bridge residue pairs at `d <= 5.0 A`, anchored by the conserved inter-subunit `ARG56-ASP144` ion pair) occur in the four Family I NTD+CTD capsomere seams (`A-A_op36`, `C-D_op1`, `B-C_op1`, `B-D_op2`), which carry **`0/10` (`0.0%`) of pi-stacking candidates**.

#### 3.3.3 Degree- and Interface-Membership-Preserving Permutation Null (`degree_and_interface_preserving_permutation_null`, `n_permutations = 10000`, `seed = 20261107`)

To account for vertex-degree heterogeneity and sub-interface contact multiplicity without assuming i.i.d. edges, `src/m7_arc_level3_interface_comparison.py` executes `10,000` deterministic permutations (`numpy.random.default_rng(20261107)`):

| Permutation Null Scheme | Preserved Graph / Interface Invariants | Observed Statistic | Permutation Null Mean ($\pm$ SD) | Fold Enrichment / $Z$-Score | Empirical Permutation $p$-Value (Phipson-Smyth `(k+1)/(N+1)`) | Epistemic Qualification |
| :--- | :--- | :---: | :---: | :---: | :---: | :---: |
| **1A. Homologous `83`-Residue Domain Permutation (Same-Domain Contacts `NTD-NTD + CTD-CTD`)** | Exact `83`-residue interfacial composition (`38` NTD, `45` CTD) and each sequence position's contact degree across `204` edges & `10` sub-interfaces | `131 / 204` same-domain (`56` NTD-NTD, `75` CTD-CTD, `73` NTD-CTD) | `103.57` ($\pm 10.90$) (`43.27` NTD-NTD, `100.43` NTD-CTD, `60.31` CTD-CTD) | $Z = +2.5157$ | **`0.006499`** (`64 / 10,000` permutations $\ge 131$; omnibus degree-adjusted $\chi^2 = 14.8162$, $p_{\text{perm}} = 0.045595$) | `:Incomplete` (exploratory screening) |
| **1B. Chain-Specific `235`-Node Domain Permutation (Same-Domain Contacts)** | Exact `235` `(chain_id, seq_id)` vertex composition (`108` NTD, `127` CTD) and per-chain vertex degrees across `204` edges | `131 / 204` same-domain | `102.23` ($\pm 7.16$) | $Z = +4.0182$ | **`0.000100`** (`0 / 10,000` permutations $\ge 131$) | `:Incomplete` (exploratory screening) |
| **1C. Full `166`-Resolved-Sequence Domain Permutation (`41..206`: `85` NTD, `81` CTD)** | Full resolved `CA` primary sequence (`85` NTD, `81` CTD) restricted to the `83` interfacial positions and their contact degrees | `131 / 204` same-domain | `103.39` ($\pm 11.02$) (`54.05` NTD-NTD, `100.61` NTD-CTD, `49.34` CTD-CTD) | $Z = +2.5059$ | **`0.008599`** (`85 / 10,000` permutations $\ge 131$) | `:Incomplete` (exploratory screening) |
| **1D. Contiguous `166`-Residue Circular-Shift Autocorrelation Null (`166` cyclic shifts $s \in \{0, \dots, 165\}$)** | Primary-sequence contiguity (`85` contiguous NTD vs `81` contiguous CTD positions) and 3D globular spatial autocorrelation (**Confound & Restraint Framing**) | `131 / 204` aggregate same-domain; `71 / 75` (`94.67%`) `CTD-CTD` across the `6` non-Family-I interfaces | `116.29` ($\pm 18.80$) aggregate same-domain | Aggregate: `38 / 166` shifts $\ge 131$ (**does NOT beat sequence contiguity**); Non-Family-I `CTD-CTD`: `3 / 166` shifts $\ge 71 / 75$ (observed `125\|126` hinge + 2 immediate neighbors $s \in \{-1, 0, +1\}$) | Aggregate: **`0.228916`** (`38 / 166`, confound); Non-Family-I `CTD-CTD` purity: **`0.018072`** (`3 / 166`, discrete floor $1/166 = 0.006024$, boundary-localization restraint) | `:Incomplete` (confound & boundary-localization restraint, not strong support) |
| **2A-i. Sub-Interface Assignment Permutation: Family II Hydrophobic Clustering** | Exact sub-interface edge counts (`55` Family II, `129` Family I, `20` Family III edges) and `20` hydrophobic contact edges | `13 / 20` edges in Family II | `5.3754` | `2.4184x` | **`0.000400`** (`3 / 10,000` permutations $\ge 13$) | `:Incomplete` (exploratory screening) |
| **2A-ii. Sub-Interface Assignment Permutation: Family II `<= 5.0 A` Pi-Stacking Clustering** | Exact sub-interface edge counts (`55` Family II edges) and `8` `d <= 5.0 A` pi-stacking contact edges (`8/10` across all cutoffs) | `8 / 8` edges in Family II | `2.1606` | `3.7027x` | **`0.000300`** (`2 / 10,000` permutations $= 8$) | `:Incomplete` (exploratory screening) |
| **2A-iii. Sub-Interface Assignment Permutation: Family I Salt-Bridge Residue-Pair Clustering** | Exact sub-interface edge counts (`129` Family I edges) and `6` salt-bridge residue pairs (`9` atom-level hits) | `6 / 6` residue pairs (`9/9` atom hits) in Family I | `3.7942` | `1.5814x` | **`0.063594`** (`635 / 10,000` permutations $= 6$) | `:Incomplete` (exploratory screening) |
| **2B-i. Homologous `83`-Residue Amino-Acid (`comp_id`) Permutation: Family II Nonpolar Pairs** | Exact `83`-residue amino-acid type multiplicity and fixed contact graph topology/degrees (`topological_graph_fixed_caveat`) | `18` nonpolar-nonpolar pairs in Family II | `6.0115` | `2.9943x` | **`0.029797`** (`297 / 10,000` permutations $\ge 18`) | `:Incomplete` (exploratory screening) |
| **2B-ii. Homologous `83`-Residue Amino-Acid (`comp_id`) Permutation: Family II Aromatic Pairs** | Exact `83`-residue aromatic multiplicity (`PHE/TYR/TRP/HIS`) and fixed contact graph topology/degrees (`topological_graph_fixed_caveat`) | `9` aromatic-aromatic pairs in Family II | `2.2859` | `3.9372x` | **`0.043496`** (`434 / 10,000` permutations $\ge 9`) | `:Incomplete` (exploratory screening) |
| **2B-iii. Homologous `83`-Residue Amino-Acid (`comp_id`) Permutation: Family I Ion Pairs** | Exact `83`-residue cationic/anionic multiplicity (`ARG/LYS` vs `ASP/GLU`) and fixed contact graph topology/degrees (`topological_graph_fixed_caveat`) | `16` opposite-charge pairs in Family I | `5.4123` | `2.9562x` | **`0.011199`** (`111 / 10,000` permutations $\ge 16`) | `:Incomplete` (exploratory screening) |

**Fixed-Graph Topological Caveat (`topological_graph_fixed_caveat`)**:
> *"This test is conditional on the fixed 204-edge contact graph and tests whether amino-acid labels are unusually distributed over already-observed edges; it does not model mutation-induced side-chain repacking or recompute the contact graph after a mutation, and its enrichment results must not be described as predictions of mutational delta-delta-G or assembly effects."*

#### 3.3.4 Explicit 7-Row Graph State-Space Specification (`state_space_specification`) and Constrained Double-Edge-Swap MCMC Audit (`constrained_graph_state_space_connectivity_and_mcmc_audit`)

Before presenting any MCMC statistic, `constrained_graph_state_space_connectivity_and_mcmc_audit` defines an explicit 7-row `state_space_specification` covering **Null Scheme I (`global_degree_preserving_simple_graph_double_edge_swap_mcmc`)** and **Null Scheme II (`within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc`)** as strictly separate null models:

| Row | `state_space_specification` Dimension | Null Scheme I: `global_degree_preserving_simple_graph_double_edge_swap_mcmc` (Global Simple Graph) | Null Scheme II: `within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc` (Within-Subinterface Bipartite) |
| :---: | :--- | :--- | :--- |
| **1** | **`vertex_identity`** | `235` chain-specific interfacial residue nodes `(chain_id in {A,B,C,D}, seq_id in 41..206, comp_id, domain_label in {NTD: 41..125, CTD: 126..206})` (`108` NTD nodes, `127` CTD nodes; `83` homologous sequence positions on the quotient graph: `38` NTD, `45` CTD). | Same `235` chain-specific interfacial residue nodes partitioned across the `10` symmetry-inequivalent sub-interfaces with fixed `(chain_id, seq_id)` and NTD/CTD domain labels. |
| **2** | **`graph_class`** | Undirected simple graph (`235` vertices, `204` edges at $d \le 5.0\text{ \AA}$, zero self-loops, zero multi-edges). | Product of `10` simple sub-interface graphs (`204` total edges, zero self-loops, zero multi-edges): within each of the `8` heterotypic sub-interfaces (`C-D_op1`, `B-C_op1`, `B-D_op2`, `A-D_op1`, `B-C_op33`, `B-D_op33`, `A-D_op36`, `A-B_op2`), the sub-graph is strictly bipartite between the two chain classes, whereas `A-A_op36` and `C-C_op33` are homotypic symmetry-interface graphs. |
| **3** | **`invariants`** | Preserves exact per-node degree sequence $\deg(v)$ across all `235` nodes, total edge count (`204`), and fixed node NTD/CTD domain labels ($S_N = 185$ NTD stubs, $S_C = 223$ CTD stubs in empirical `6TAP` at `NTD: 41..125, CTD: 126..206`; or $S_N = 169$ NTD stubs, $S_C = 239$ CTD stubs in the illustrative `35`-stub-difference parameterization), while allowing edges to rewire across sub-interfaces. | Preserves exact sub-interface membership (`10` sub-interfaces with empirical `6TAP` edge counts `41, 34, 28, 28, 27, 26, 7, 7, 5, 1`; illustrative partition `41, 34, 34, 20, 21, 34, 7, 6, 6, 1`), exact within-subinterface per-node degree sequence $\deg_I(v)$, and chain-bipartite partition on the `8` heterotypic interfaces. |
| **4** | **`swap_constraints`** | Propose two distinct edges $e_1 = (u, v), e_2 = (x, y)$ and candidate rewiring $e_1' = (u, y), e_2' = (x, v)$ or $(u, x), (v, y)$. Reject if $e_1'$ or $e_2'$ creates a self-loop ($u=y$ or $x=v$) or duplicates an existing edge (multi-edge collision). **Detailed balance / self-loop persistence**: on rejection, retain $G_{t+1} = G_t$ (lazy Metropolis-Hastings persistence, acceptance rate `0.9796`), ensuring symmetric transition probabilities $P(G \to G') = P(G' \to G)$. | Propose two distinct edges $e_1 = (u, v), e_2 = (x, y)$ within the same sub-interface $I$ (`2,865` unordered within-subinterface pairs) and candidate chain-preserving rewiring $e_1' = (u, y), e_2' = (x, v)$. Reject if self-loop or multi-edge collision occurs. **Detailed balance / self-loop persistence**: on rejection, retain $G_{t+1} = G_t$ (acceptance rate `0.8160`), ensuring $P(G \to G') = P(G' \to G)$. |
| **5** | **`stationary_distribution`** | Uniform distribution $\pi(G) = \text{const}$ over the connected component of reachable `235`-vertex, `204`-edge simple graphs with exact degree sequence $\deg(v)$ (verified via symmetric proposal + persistence on rejection). | Uniform distribution $\pi(G) = \text{const}$ over the connected component of reachable simple sub-interface graphs with exact within-subinterface degree sequences $\deg_I(v)$ and chain-bipartite partitions. |
| **6** | **`connectivity`** | Connectedness of simple graphs with a given degree sequence under double-edge swaps (Taylor 1981; Ryser 1957; Fosdick et al. 2018 *SIAM Rev.*), qualified by whether additional chain/spatial restrictions are imposed ($E_{\text{same}}$ explores `[77, 135]`, mean `102.6735 +/- 7.0909`, $z = +3.9948$, `46 / 200,000` steps $\ge 131$, $p_{\text{PS}} = 0.000235$). In $G_0$, out of $\binom{204}{2} = 20,706$ unordered edge pairs (`41,412` directed proposals), **`20,441 / 20,706` (`98.72%`)** unordered pairs admit $\ge 1$ valid simple swap (**`40,616 / 41,412` (`98.08%`)** of directed proposals), and **`6,753 / 20,706` (`32.61%` of all unordered pairs; `6,753 / 20,441 = 33.04%` of valid unordered pairs; `10,953 / 41,412 = 26.45%` of all directed proposals; `10,953 / 40,616 = 26.97%` of valid directed proposals)** alter $E_{\text{same}}$ by $\pm 2$. | **Structural Confound & Restraint (NOT Positive Support)**: **6 of the 10 sub-interfaces** (`A-D_op1`, `B-C_op33`, `C-C_op33`, `B-D_op33`, `A-D_op36`, `A-B_op2`, totaling `71 / 75` pure `CTD-CTD` contacts) have either `0` NTD stubs (`4` sub-interfaces: `C-C_op33`, `B-D_op33`, `A-D_op36`, `A-B_op2`, and also `A-D_op1, C-C_op33, B-D_op33, A-B_op2` when boundary `ALA125` is grouped with the hinge) or single-boundary `ALA125` NTD stubs that cannot reduce `CTD-CTD` below `71/75`; in addition, `B-D_op2` has `28 NTD / 0 CTD` stubs on Chain D (`0` color-changing swaps). Out of `2,865` unordered within-subinterface edge pairs (`2,377 / 2,865 = 82.97%` valid swaps), only `48` swaps (**`48 / 2,865 = 1.68%` of all unordered within-subinterface pairs; `48 / 2,377 = 2.02%` of valid within-subinterface swaps**) alter domain-pair counts, trapping $E_{\text{same}} \in [123, 135]$ (mean `128.4536 +/- 1.8941`, $p_{\text{PS}} = 0.209264$) and proving that once sub-interface membership is conditioned on, domain identity adds almost no signal beyond the sub-interface boundary itself. |
| **7** | **`mixing`** | `4` independent chains (`seed = 20261108 + c`), `5,000` burn-in, `50,000` post-burn-in steps/chain (`200,000` total steps): $\rho_1 = 0.9903$, $\rho_{10} = 0.9066$, $\rho_{50} = 0.6117$, $\tau_{\text{int}} = 209.61$, $\text{ESS}_{\text{total}} = 954.15$ (`238.54`/chain, `954` effectively independent states), Gelman-Rubin $\hat{R} = 1.00259$, $\text{MCSE}_{\text{mean}} = 0.2296$, $p_{\min,\text{ESS}} = 1/(\text{ESS}+1) = 0.001047$. | `4` independent chains (`seed = 20261108 + c`), `5,000` burn-in, `50,000` post-burn-in steps/chain (`200,000` total steps): $\rho_1 = 0.9891$, $\rho_{10} = 0.8962$, $\rho_{50} = 0.5754$, $\tau_{\text{int}} = 193.95$, $\text{ESS}_{\text{total}} = 1,031.18$ (`257.79`/chain, `1,031` effectively independent states), Gelman-Rubin $\hat{R} = 1.00159$, $\text{MCSE}_{\text{mean}} = 0.0590$, $p_{\min,\text{ESS}} = 1/(\text{ESS}+1) = 0.000969$. |

##### Per-Subinterface Swap Mobility Under Null Scheme II (`per_subinterface_swap_mobility`: Both Denominators Explicitly Reported)

| Sub-Interface Class | Family | Edges ($5\text{ \AA}$) | Chain $i$ Stubs (`NTD / CTD`) | Chain $j$ Stubs (`NTD / CTD`) | Total Stubs (`NTD / CTD`) | Observed (`NN / NC / CC`, $E_{\text{same}}$) | Unordered Edge Pairs $\binom{m_I}{2}$ | Valid Simple Swaps (`% of Unordered Pairs`) | Color-Changing Swaps (`% of Unordered Pairs` / `% of Valid Swaps`) | Theoretical $E_{\text{same}}$ Corridor | Frozen Under Scheme II? |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `A-A_C5_pentamer_op36` | I | `41` | `20 / 21` | `40 / 1` | `60 / 22` | `20 / 20 / 1` (`21`) | `820` | `695` (`84.76%`) | `20` (`2.44%` / `2.88%`) | `[19, 21]` | No |
| `C-D_C1_intra_asu_op1` | I | `34` | `33 / 1` | `12 / 22` | `45 / 23` | `12 / 21 / 1` (`13`) | `561` | `490` (`87.34%`) | `12` (`2.14%` / `2.45%`) | `[11, 13]` | No |
| `B-C_C3_inter_asu_hexamer_op33` | II | `28` | `1 / 27` | `1 / 27` | `2 / 54` (`0 / 56` excl. `ALA125`) | `0 / 2 / 26` (`26`) | `378` | `303` (`80.16%`) | `1` (`0.26%` / `0.33%`) | `[26, 28]` | At min `CTD-CTD = 26` |
| `B-D_C2_dimer_op2` | I | `28` | `17 / 11` | `28 / 0` | `45 / 11` | `17 / 11 / 0` (`17`) | `378` | `324` (`85.71%`) | `0` (`0.00%` / `0.00%`) | `[17, 17]` | **Yes (`0` color swaps)** |
| `A-D_C1_intra_asu_op1` | II | `27` | `1 / 26` | `1 / 26` | `2 / 52` (`0 / 54` excl. `ALA125`) | `0 / 2 / 25` (`25`) | `351` | `269` (`76.64%`) | `1` (`0.28%` / `0.37%`) | `[25, 27]` | At min `CTD-CTD = 25` |
| `B-C_C1_intra_asu_op1` | I | `26` | `24 / 2` | `7 / 19` | `31 / 21` | `7 / 17 / 2` (`9`) | `325` | `279` (`85.85%`) | `14` (`4.31%` / `5.02%`) | `[5, 9]` | No |
| `C-C_C3_trimer_op33` | III | `7` | `0 / 7` | `0 / 7` | `0 / 14` | `0 / 0 / 7` (`7`) | `21` | `10` (`47.62%`) | `0` (`0.00%` / `0.00%`) | `[7, 7]` | **Yes (`0` NTD stubs)** |
| `B-D_C3_inter_asu_vertex_op33` | III | `7` | `0 / 7` | `0 / 7` | `0 / 14` | `0 / 0 / 7` (`7`) | `21` | `6` (`28.57%`) | `0` (`0.00%` / `0.00%`) | `[7, 7]` | **Yes (`0` NTD stubs)** |
| `A-D_C5_inter_asu_seam_op36` | III | `5` | `0 / 5` | `0 / 5` | `0 / 10` | `0 / 0 / 5` (`5`) | `10` | `1` (`10.00%`) | `0` (`0.00%` / `0.00%`) | `[5, 5]` | **Yes (`0` NTD stubs)** |
| `A-B_C2_peripheral_op2` | III | `1` | `0 / 1` | `0 / 1` | `0 / 2` | `0 / 0 / 1` (`1`) | `0` | `0` (`0.00%`) | `0` (`0.00%` / `0.00%`) | `[1, 1]` | **Yes (`0` NTD stubs)** |
| **Total Across 10 Sub-Interfaces** | **All** | **`204`** | **`96 / 108`** | **`89 / 115`** | **`185 / 223`** | **`56 / 73 / 75` (`131`)** | **`2,865`** | **`2,377` (`82.97%`)** | **`48` (`1.68%` of pairs / `2.02%` of valid)** | **`[123, 135]`** | **5 frozen; 6 non-Family-I $\ge 71/75$ `CTD-CTD`** |

> **Reconciliation of Earlier Stale Prose Table (`stale_table_reconciliation_note`)**: An earlier prose summary mixed the illustrative partition's edge/pair counts (`41, 34, 34, 20, 21, 34, 7, 6, 6, 1` $\to$ `190` pairs for `B-D_C2`, `210` pairs for `A-D_C1`, summing to `2,954`) with the empirical `6TAP` valid swap counts (`324` from $m=28, \binom{28}{2}=378$ pairs for `B-D_C2`, and `269` from $m=27, \binom{27}{2}=351$ pairs for `A-D_C1`, summing to `2,377` valid swaps out of `2,865` empirical unordered pairs). In empirical `6TAP`, every sub-interface strictly satisfies `valid_simple_bipartite_swaps <= unordered_edge_pairs` (`324 <= 378` and `269 <= 351`). Likewise, reporting both denominators (`48 / 2,865 = 1.68%` of unordered within-subinterface pairs and `48 / 2,377 = 2.02%` of valid within-subinterface swaps, alongside `6,753 / 20,706 = 32.61%` of unordered global pairs and `6,753 / 20,441 = 33.04%` of valid global swaps [`10,953 / 41,412 = 26.45%` of directed proposals, `10,953 / 40,616 = 26.97%` of valid directed proposals]) keeps all swap-mobility comparisons in identical units.

##### Precise Qualification of the Stub-Difference Identity (`stub_difference_conservation_identity`)

- **Role**: `"CONSERVATION_INVARIANT_AND_IMPLEMENTATION_CHECK_NOT_EMPIRICAL_SIGNIFICANCE_STATISTIC"`
- **Algebraic Proof**:
  Let $S_C = 2 E_{CC} + E_{CN}$ denote the total number of CTD half-edge stubs and $S_N = 2 E_{NN} + E_{CN}$ denote the total number of NTD half-edge stubs across the `204` edges ($S_C + S_N = 408$). Subtracting $S_N$ from $S_C$ eliminates $E_{CN}$ identically in every degree-and-color-preserving graph realization:
  $$S_C = 2 E_{CC} + E_{CN}, \quad S_N = 2 E_{NN} + E_{CN} \implies E_{CC} - E_{NN} = \frac{S_C - S_N}{2}$$
  Under the pre-run illustrative parameterization ($S_C = 239, S_N = 169$, corresponding to $E_{CC} = 83, E_{NN} = 48, E_{CN} = 73$ at $E_{\text{same}} = 131$ where `35` equals the count of 4-chain conserved sequence positions):
  $$S_C = 2 E_{CC} + E_{CN} = 239, \quad S_N = 2 E_{NN} + E_{CN} = 169 \implies E_{CC} - E_{NN} = \frac{S_C - S_N}{2} = \frac{239 - 169}{2} = 35$$
  Under the empirical `6TAP` $d \le 5.0\text{ \AA}$ partition (`NTD: 41..125`, `CTD: 126..206`, where $E_{CC} = 75, E_{NN} = 56, E_{CN} = 73$):
  $$S_C = 2(75) + 73 = 223, \quad S_N = 2(56) + 73 = 185 \implies E_{CC} - E_{NN} = \frac{223 - 185}{2} = \frac{38}{2} = 19$$
  In both cases, $E_{CC} - E_{NN}$ is an exact mathematical consequence of the degree-stub difference (`239 - 169 = 70` giving `35`, or `223 - 185 = 38` giving `19`, verified with `min = 19, max = 19, variance = 0.0` across all `200,000` MCMC steps in both Scheme I and Scheme II) rather than an empirical significance statistic, whereas $E_{\mathrm{same}} = E_{CC} + E_{NN} = 204 - E_{CN}$ varies in steps of $\pm 2$.

##### Distinct Roles of the Methodological References & Monte Carlo / MCMC P-Value Resolution

1. **Polaris (`arXiv:2409.01363`)**: Defines the colored multigraph configuration model preserving the joint color matrix and vertex degree sequence; transferring this framework to a simple, chain-bipartite protein contact graph requires simple-graph collision rejection (zero self-loops, zero multi-edges) and chain-bipartite endpoint constraints, and when the $2 \times 2$ domain color matrix is relaxed to test domain assortativity, marginal stub counts $S_C, S_N$ still lock $E_{CC} - E_{NN} = (S_C - S_N)/2$.
2. **Clauset, Fosdick, & Dutta (`arXiv:2105.12120` / Fosdick et al. 2018 *SIAM Rev.*)**: Establishes stub-labeled vs vertex-labeled simple-graph state spaces, proves that lazy self-loop persistence ($G_{t+1} = G_t$ on rejected swaps) is required for a uniform stationary distribution, and defines the MCMC convergence diagnostics ($\rho_k, \tau_{\mathrm{int}}, \mathrm{ESS}, \hat{R}$) used here.
3. **Topologically Invariant Permutation Test (`arXiv:2511.06153`)**: Persistence-diagram / Wasserstein permutation testing on functional brain networks; cited solely as a domain comparator for topology-invariant permutation testing and **explicitly NOT a proof of irreducibility or mixing for the M7 double-edge-swap kernel** (`is_proof_of_m7_double_edge_swap_irreducibility_or_mixing = false`).
4. **KDD 2026 (*Testing Statistical Dependence in Labeled Graphs under Mismatches*)** and **MANGO (2026 preprint)**: Formulate topology-and-label dependence testing under graph mismatches and degree-aware clustering comparators conditioning on fixed topology and degree composition.
5. **Phipson & Smyth (2010, `doi:10.2202/1544-6115.1585`) & Effective Sample Floor**: For i.i.d. permutations ($B = 10,000$), the minimum attainable resolution floor is $p_{\min,\mathrm{perm}} = \frac{1}{B + 1} = \frac{1}{10001} \approx 9.999 \times 10^{-5}$. For correlated MCMC chains (`200,000` raw post-burn-in steps across 4 chains), `monte_carlo_and_mcmc_p_value_resolution` reports both the raw step count and the **effective independent sample floor** $p_{\min,\mathrm{ESS}} = \frac{1}{\mathrm{ESS} + 1}$ (`0.001047` for Scheme I at $\mathrm{ESS} = 954.15$; `0.000969` for Scheme II at $\mathrm{ESS} = 1,031.18$), alongside $\mathrm{MCSE}$ and exact 95% Clopper-Pearson and Wilson score intervals on $\hat{p}$ (for Scheme I `46 / 200,000`: $\hat{p}_{\mathrm{PS}} = 0.000235$, $\mathrm{MCSE}_{\hat{p}} = 0.000496$, Clopper-Pearson 95% CI `[0.000168, 0.000307]`, Wilson 95% CI `[0.000172, 0.000307]`; for Scheme II `41,852 / 200,000`: $\hat{p}_{\mathrm{PS}} = 0.209264$, $\mathrm{MCSE}_{\hat{p}} = 0.012665$, Clopper-Pearson 95% CI `[0.207479, 0.211049]`, Wilson 95% CI `[0.207483, 0.211048]`).

#### 3.3.5 Multiplicity Adjustment Analysis Across Hypothesis Families A and B (`multiplicity_adjustment_analysis`)

| Test ID | Hypothesis Family | Exceedances ($k / N$) | Unadjusted $p$ | Within-Family Holm-Bonferroni $p_{\mathrm{Holm,fam}}$ | Within-Family BH-FDR $q_{\mathrm{BH,fam}}$ | 12-Test Omnibus Holm-Bonferroni $p_{\mathrm{Holm,12}}$ | 12-Test Omnibus BH-FDR $q_{\mathrm{BH,12}}$ | Survives Within-Family (`Holm` / `BH` at `0.05`) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| `chain_specific_235_node_same_domain` | Family A (`m = 6`) | `0 / 10,000` | `0.000100` | **`0.000600`** | **`0.000600`** | **`0.001200`** | **`0.001200`** | **Yes / Yes** |
| `homologous_83_residue_same_domain` | Family A (`m = 6`) | `64 / 10,000` | `0.006499` | **`0.032497`** | **`0.017198`** | `0.058494` | **`0.019498`** | **Yes / Yes** |
| `full_166_resolved_sequence_same_domain` | Family A (`m = 6`) | `85 / 10,000` | `0.008599` | **`0.034397`** | **`0.017198`** | `0.068793` | **`0.020638`** | **Yes / Yes** |
| `circular_shift_166_non_family_I_ctd_purity` | Family A (`m = 6`) | `3 / 166` | `0.018072` | `0.054217` | **`0.027108`** | `0.108434` | **`0.030981`** | No / **Yes** |
| `homologous_83_residue_omnibus_chi_square` | Family A (`m = 6`) | `455 / 10,000` | `0.045595` | `0.091191` | `0.054715` | `0.173983` | `0.054715` | No / No |
| `circular_shift_166_aggregate_same_domain` | Family A (`m = 6`) | `38 / 166` | `0.228916` | `0.228916` | `0.228916` | `0.228916` | `0.228916` | No / No |
| `subinterface_edge_permutation_family_II_pi_stacking` | Family B (`m = 6`) | `2 / 10,000` | `0.000300` | **`0.001800`** | **`0.001200`** | **`0.003300`** | **`0.001600`** | **Yes / Yes** |
| `subinterface_edge_permutation_family_II_hydrophobic` | Family B (`m = 6`) | `3 / 10,000` | `0.000400` | **`0.002000`** | **`0.001200`** | **`0.004000`** | **`0.001600`** | **Yes / Yes** |
| `homologous_83_comp_id_permutation_family_I_ion_pair` | Family B (`m = 6`) | `111 / 10,000` | `0.011199` | **`0.044796`** | **`0.022398`** | `0.078392` | **`0.022398`** | **Yes / Yes** |
| `homologous_83_comp_id_permutation_family_II_nonpolar` | Family B (`m = 6`) | `297 / 10,000` | `0.029797` | `0.089391` | **`0.044696`** | `0.148985` | **`0.044696`** | No / **Yes** |
| `homologous_83_comp_id_permutation_family_II_aromatic` | Family B (`m = 6`) | `434 / 10,000` | `0.043496` | `0.089391` | `0.052195` | `0.173983` | `0.054715` | No / No |
| `subinterface_edge_permutation_family_I_salt_bridge` | Family B (`m = 6`) | `635 / 10,000` | `0.063594` | `0.089391` | `0.063594` | `0.173983` | `0.069375` | No / No |

---

## 4. Three Physical RNA Arrangements & Four-Branch Literature Separation

### 4.1 Explicit Separation of Three Distinct Physical RNA Arrangements

To prevent conflation between native capsid structure, engineered delivery vehicles, and non-encapsulated extracellular vesicles, Level 3 explicitly separates three physical RNA arrangements:

| Arrangement ID | Physical RNA Topology | Primary Reference | Deposited mmCIF / Molecular Architecture | Epistemic Status |
| :--- | :--- | :--- | :--- | :---: |
| **`arrangement_1_internal_nc_zn_rna_binding_architecture_in_6tap`** | **Internal `NC/Zn` RNA-binding architecture in `6TAP`** | Erlendsson et al. (*Nat. Neurosci.* 2020, `PDB 6TAP`) | Polymer `_entity.id = 1` (`dArc1` polypeptide) strand `label_asym_id = E` (`auth_asym_id = E`, residues `224..250`, `NC_ZINC_KNUCKLE_FRAGMENT`) coordinating non-polymer ligand `_entity.id = 2` (`ZN`, `label_asym_id = F`, `auth_asym_id = E`, `_atom_site.auth_seq_id = 301`, `pdb_seq_num = 301`, `auth_seq_num = 251`, `ZN_COORDINATION_LIGAND`) projecting inward beneath the `C2` `B-D` hexamer seam (`r = 163.70 A` vs `CA` chains `A..D` at `171.12..180.38 A`, `E-E_op2` dimer with `12` unordered residue pairs / `21` directed contacts at `<= 5.0 A`, `Zn-Zn = 14.8186 A`). | `:Accepted` |
| **`arrangement_2_direct_lumenal_rna_cargo_encapsulation_engineered_darc1`** | **Direct lumenal RNA cargo encapsulation inside engineered `dArc1` capsids** | Lash et al. (*Nature Communications*, Nov 4, 2026, `doi:10.1038/s41467-026-79058-8`) | Reconstituted engineered `dArc1` capsids fused to heterologous RNA-binding domains (`RBDs`) that package exogenous mRNA/editor cargo directly inside the closed icosahedral capsid lumen for `SORCS2`-mediated mammalian delivery (`<= 18%` Cas9 exon skipping in `mdx` mice). | `:Accepted` (engineered system) / `:Incomplete` (native `6TAP` sub-interface release mechanics) |
| **`arrangement_3_non_encapsulated_arc_mbnl1_ev_cotransport`** | **Non-encapsulated `Arc-Mbnl1` extracellular-vesicle (`EV`) co-transport** | Zinter et al. (*bioRxiv* / `PMC13061017`, Mar 10, 2026) | Activity-dependent cellular mammalian `Arc-Mbnl1` mRNA association in membrane-enveloped extracellular vesicles where `Mbnl1` mRNA is **not** directly encapsulated inside a closed icosahedral `Arc` capsid shell. | `:Incomplete` |

### 4.2 Four-Branch Literature & Provenance Separation Table

| Branch ID | Literature / Evidence Source | Epistemic Status | Supported Claims (Within Scope) | Quarantined / Unsupported Extrapolations |
| :--- | :--- | :---: | :--- | :--- |
| **Branch 1: Native `dArc1` Structure (`6TAP`)** | Erlendsson et al. (*Nat. Neurosci.* 2020, `PDB 6TAP`); Cottee et al. (*Nat. Neurosci.* 2020); Hallin et al. (2021) | `:Accepted` | `T=4` icosahedral capsid architecture in icosahedral point group `I` (`_em_single_particle_entity.point_symmetry = "I"`, biological assembly `1`, `240` CA subunits, `12` pentamers, `30` hexamers); `204` unique residue contacts at `5.0 A` across `10` sub-interfaces (`9,570.22 A^2` representative two-sided BSA and `4,785.10 A^2` one-side area per ASU; `574,213.20 A^2` whole-capsid two-sided BSA and `287,106.60 A^2` whole-capsid one-side interface area); exact mmCIF entity reconciliation (`_entity.id = 1` chains `A..E`, `_entity.id = 2` `ZN` chain `F` with `auth_asym_id = E`). | Does **not** establish the kinetic assembly pathway, per-interface $\Delta\Delta G$, or in-vivo synaptic signaling role of individual sub-interfaces without mutational perturbation. |
| **Branch 2: 2026 Engineered `dArc1` Delivery** | Lash et al., *Nature Communications* (published Nov 4, 2026, `doi:10.1038/s41467-026-79058-8`) | `:Accepted` (engineered system) / `:Incomplete` (native `6TAP` interface function) | In-vitro reconstitution of engineered `dArc1` capsids; direct lumenal RNA cargo encapsulation via engineered RNA-binding domains (`RBD`); identification of mammalian `SORCS2` cell-surface receptor; up to `18%` Cas9 exon skipping in `mdx` mice. | Engineered `RBD` cargo tethering and mammalian `SORCS2` receptor usage do **not** identify which native `6TAP` CA sub-interface mediates physiological disassembly or native `dArc1` mRNA transfer in *Drosophila* (`:Incomplete`). |
| **Branch 3: 2026 `Arc-Mbnl1` RNA Association Preprint** | Zinter et al., *bioRxiv* / `PMC13061017` (posted Mar 10, 2026) | `:Incomplete` | Reports activity-dependent mammalian `Arc-Mbnl1` mRNA association and possible extracellular-vesicle (`EV`) co-transport where `Mbnl1` is **not** directly encapsulated inside an Arc capsid shell. | Unrefereed preprint (`:Incomplete`); explicitly distinguishes non-encapsulated membrane-enveloped `EV` co-transport from direct lumenal encapsulation inside the `6TAP` icosahedral capsid shell. |
| **Branch 4: Candidate Mutational Hotspots (Pre-Registered Hypotheses)** | M7 Level-2/Level-3 `6TAP` Residue-Interface Census (`2026-11-07`) | `:Incomplete` | Nominates 4 structurally specific residue motifs (`PHE133`, `ARG49/ASP66/GLU76/ARG111 + ARG56/ASP144`, `HIS182/HIS186`, `CYS227/CYS230/HIS235/CYS240`) classified by `perturbation_category` and `required_experimental_controls` with pre-stated experimental falsification readouts. | Strictly pre-registered structural hypotheses (`:Incomplete`); zero thermodynamic hotspot claims are promoted to `:Accepted` prior to wet-lab SEC-MALS, DSF, cryo-EM, ICP-MS/PAR, CD/NMR, and RNA-binding assays. |

### 4.3 Pre-Registered Candidate Mutational Hotspots, Perturbation Categories, Controls & Falsification Readouts (Branch 4)

| Hotspot ID | Nominated Residue Motif | Target Sub-Interface(s) | `perturbation_category` | `required_experimental_controls` | Structural Basis & Pre-Registered Falsification Readout (`:Incomplete`) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| `hotspot_1_PHE133_ctd_aromatic_lock` | `PHE133` (with `TYR126 / TYR129 / HIS162 / PHE172` aromatic cluster) | `A-D_C1_intra_asu_op1`, `B-C_C3_inter_asu_hexamer_op33` | `CA_CA_SUBUNIT_INTERFACE_ENERGETICS_PERTURBATION` | Monomer secondary/tertiary fold verification (far-UV CD and 1H-15N HSQC NMR of isolated CTD), SEC-MALS oligomerization stoichiometry, DSF thermal unfolding (`Tm`), and negative-stain/cryo-EM `T=4` capsid shell morphology. | `A:PHE133-D:PHE133` (`2.8699 A` min heavy-atom distance, `4.217 A` ring centroid) and `B:PHE133-C:PHE133` (`3.5768 A` min distance, `4.933 A` centroid) form the central inter-subunit aromatic lock across both Family II CTD dimers (`8/10` of all capsid pi-stacking candidates). **Mutation**: `F133A` / `F133E`. **Falsified if**: retains wild-type `T=4` capsid yield and hydrodynamic radius by SEC-MALS / negative-stain EM ($\Delta T_m < 1.5\text{ K}$ by DSF). |
| `hotspot_2_ARG49_ASP66_GLU76_ARG111_pentamer_salt_bridge_network` | `ARG49 / ASP66 / GLU76 / ARG111` (with conserved `ARG56-ASP144` and `GLU117-ARG199`) | `A-A_C5_pentamer_op36` (and `C-D_op1`, `B-C_op1`, `B-D_op2`) | `CA_CA_SUBUNIT_INTERFACE_ENERGETICS_PERTURBATION` | Monomer fold verification by CD/DSF across `50-500 mM NaCl` ionic-strength titration, SEC-MALS pentamer/hexamer ratio quantification, and compensatory double-mutant (`R56E/D144R`) epistasis rescue controls. | `A-A_C5_pentamer_op36` (`1,882.00 A^2`, `41` pairs) in Family I pairs conserved inter-subunit salt bridges `ARG56-ASP144` and `ARG199-GLU117` with the NTD stalk/bundle charged network. **Mutation**: `R56E/D144R`, `R49E/E76R`, `D66R/R111D` + `500 mM NaCl` titration. **Falsified if**: charge swap or high-salt screening fails to perturb pentamer/hexamer seam closure by cryo-EM/SEC-MALS, or double charge-swap fails to rescue single mutants. |
| `hotspot_3_HIS182_HIS186_ctd_trimer_contact` | `HIS182 / HIS186` (with `LYS181 / GLN197 / ASN204`) | `C-C_C3_trimer_op33`, `B-D_C3_inter_asu_vertex_op33` | `CA_CA_SUBUNIT_INTERFACE_ENERGETICS_PERTURBATION` | pH-dependent monomer stability controls (CD/DSF across `pH 5.5-8.0`), SEC-MALS capsid mass distribution, and RNase protection / capsid permeability assays. | `C:HIS182-C:HIS182` (`6.893 A` centroid) and `B:HIS182-D:HIS182` (`6.808 A` centroid) form symmetric `C3` vertex contacts in the `Family_III_Residual_Distal_CTD_Hierarchy` Tier 1 `C3` core (`Jaccard = 0.5000`, `Sorensen-Dice = 0.6667`). **Mutation**: `H182A/H186A` + `pH 5.5..8.0` titration. **Falsified if**: capsid porosity, RNase protection, and `C3` vertex stability are invariant across `pH 5.5..8.0` and unaffected by `H182A/H186A`. |
| `hotspot_4_CYS227_CYS230_HIS235_CYS240_nc_zn_knuckle` | `CYS227 / CYS230 / HIS235 / CYS240` (`CCHC` Zn-knuckle in `NC` chain `E`, `224..250`) | `E-E_C2_nc_dimer_op2`, `NC_chain_E_ZN_F301_coordination_sphere` | `METAL_COORDINATION_AND_NC_FOLDING_PERTURBATION_NOT_ORDINARY_CA_INTERFACE_HOTSPOT` | **Zn occupancy (`ICP-MS/PAR`), `NC` folding (`CD/NMR`), and RNA-binding controls (`EMSA`/fluorescence polarization with `dArc1` 3'-UTR RNA) alongside empty-capsid `CA` shell assembly controls (`SEC-MALS`/cryo-EM)** so `NC` misfolding or $\text{Zn}^{2+}$ loss is never confounded with ordinary `CA-CA` interface energetics. | Interior `NC` strand `E` (`label_asym_id = E`, `label_entity_id = "1"`, `224..250`) coordinates `ZN` (`label_asym_id = F`, `auth_asym_id = E`, `label_entity_id = "2"`, `_atom_site.auth_seq_id = 301`, `pdb_seq_num = 301`, `auth_seq_num = 251`) with tetrahedral geometry (`tau_4 = 0.9331`, mean angle `109.45 deg`) and forms a `C2` dimer across `op2` (`12` unordered residue pairs / `21` directed contacts at `<= 5.0 A`, `Zn-Zn` distance `14.8186 A`). **Mutation**: `C227S/C230S` or `TPEN/EDTA` $\text{Zn}^{2+}$ depletion. **Falsified if**: verified $\text{Zn}^{2+}$ depletion causes no change in `dArc1` 3'-UTR binding stoichiometry or `C2` `NC`-dimer ordering while preserving empty-capsid `CA` shell assembly. |

---

## 5. Reconciliation and Restoration of the `2,996 -> 3,059` (`+63` Individual `@test` Assertions) Julia Test Count (`julia_assertion_count_reconciliation_3059_to_2996` and `julia_assertion_count_restoration_2996_to_3059`) and Complete 7-Number Assertion Audit Trail

To maintain complete audit transparency across `VERIFICATION_LEDGER_2026-11-02.json`, `docs/VERIFICATION_LEDGER_2026-11-02.md`, and this specification, both the historical `-63` consolidation (`3,059 -> 2,996`) and the Commit 5 `+63` restoration (`2,996 -> 3,059`) of individual `@test` statements are recorded below, followed by the complete 7-number cross-repository assertion audit trail (`2,996`, `3,059`, `3,830`, `3,893`, `3,913`, `4,017`, `4,080`):

| Julia Test Suite File | Consolidated Baseline (`2,996`) | Restored Commit 5 Runtime `@test` Total (`3,059`) | Restored File-Level Static `@test` Macro Count | Net Runtime Assertion Delta | Restoration & Reconciliation Explanation |
| :--- | :---: | :---: | :---: | :---: | :--- |
| `code/test_FlowCytometryObservationModel.jl` | `293` | `293` | `285` | `0` | Unchanged (`293 / 293` runtime `@test` passes across `285` static `@test` macros). |
| `code/test_MicrofluidicClassificationEngine.jl` | `1,682` | `1,735` | `764` (`711 + 53`) | `+53` | Unrolled `53` `L89..L103` predicates previously consolidated via `&&` in `l89_103_ok` back into `53` individual `@test` / `@test_throws` statements (`711 -> 764` static macros, `1,682 -> 1,735` runtime passes) so any individual failure is immediately localized by line number. |
| `code/test_NuclearReactiveTransportBridge.jl` | `332` | `332` | `276` | `0` | Unchanged (`332 / 332` runtime `@test` passes across `276` static `@test` macros). |
| `code/test_M7VisualizationBridge.jl` | `689` | `699` | `582` (`236` in v1.0 / ALARA `@testset` block) | `+10` | Unrolled `10` `L89..L103` predicates previously consolidated via `&&` back into `10` individual `@test` statements (`226 -> 236` static `@test` macros in the v1.0 / ALARA `@testset` block, `572 -> 582` file-level static macros, `689 -> 699` runtime passes). |
| **Repository Julia Total** | **`2,996`** | **`3,059`** | **`1,907`** (`236` in v1.0 / ALARA block) | **`+63`** | **All `3,059 / 3,059` (`293 + 1,735 + 332 + 699`) individual runtime assertions execute and pass (`3,887` across all Julia suites, `3,893` with Python scene adapter, `4,080` grand total with NeuralCompose).** |

### 5.1 Complete 7-Number Cross-Repository Assertion Audit Trail (`2,996`, `3,059`, `3,830`, `3,893`, `3,913`, `4,017`, `4,080`)

| Reported Assertion Total | Scope & Suite Composition | Arithmetic Breakdown | Why This Number Appeared in Earlier Snapshots |
| :---: | :--- | :--- | :--- |
| **`2,996`** | `flow-cytometry-model-bridge` (4 Julia suites, consolidated) | `293 + 1,682 + 332 + 689 = 2,996` | Intermediate snapshot where `63` `L89..L103` predicates (`53` in `test_MicrofluidicClassificationEngine.jl`, `10` in `test_M7VisualizationBridge.jl`) were conjugated via `&&` into composite `@test` macros. |
| **`3,059`** | `flow-cytometry-model-bridge` (4 Julia suites, restored) | `293 + 1,735 + 332 + 699 = 3,059` (`2,996 + 63`) | Current restored total with all `63` `L89..L103` predicates unrolled into individual top-level `@test` / `@test_throws` statements across `1,907` static `@test` macros. |
| **`3,830`** | 7-Suite Cross-Repo Total (under consolidated `2,996`) | `2,996 (bridge Julia) + 327 (M7Multiscale) + 501 (SMBR.jl) + 6 (Python Sobolev) = 3,830` | 7-suite cross-repository total before unrolling the `63` conjugated `&&` predicates. |
| **`3,893`** | 7-Suite Cross-Repo Total (under restored `3,059`) | `3,059 (bridge Julia) + 327 (M7Multiscale) + 501 (SMBR.jl) + 6 (Python Sobolev) = 3,893` (`3,830 + 63`) | Current 7-suite cross-repository total (`3,887` Julia + `6` Python Sobolev adapter). |
| **`3,913`** | 7-Suite Cross-Repo Total + 20 Pre-MCMC ARC Python Tests | `3,893 (7-suite total) + 20 (12 Level-2 + 8 pre-MCMC Level-3 pytest) = 3,913` | Intermediate total before adding the `2` permutation/MCMC Level-3 Python `pytest` checks (with all `22` current ARC Python tests + `6` Sobolev Python tests, the sum is `3,893 + 22 = 3,915`). |
| **`4,017`** | 8-Suite Cross-Repo Grand Total (under consolidated `2,996`) | `3,830 (7-suite consolidated) + 187 (NeuralCompose Tests/eval) = 4,017` | 8-suite grand total including `NeuralCompose` `Tests/eval` (`187` tests across 22 files) under the consolidated `2,996` bridge Julia count. |
| **`4,080`** | 8-Suite Cross-Repo Grand Total (under restored `3,059`) | `3,893 (7-suite restored) + 187 (NeuralCompose Tests/eval) = 4,080` (`4,017 + 63`) | Current restored 8-suite cross-repository grand total (excluding the `22` M7 ARC Level-2/3 Python `pytest` checks; `4,080 + 22 = 4,102` when included). |


