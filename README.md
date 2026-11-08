# ARC Capsid Geometry & AHRT v0.1 Multi-Platform Visualization Repository

**Epistemic Status**: `[MEASURED GEOMETRY & TOPAS-NBIO TRACK SERIES - NO MECHANISM CLAIMED]`

This repository unifies the structural geometry, symmetry-normalized residue interface graphs (M7 ARC Levels 1 to 3), and **Arc-Hypothesis Radiation-Transport (`AHRT v0.1`)** radiolysis workbench into a single portable package for **UCSF ChimeraX 1.12**, **PyMOL**, **MATLAB R2026b**, and **RStudio / R**.

---

## 1. Deposited Structural Coordinates (`structures/`)

All structural coordinates in `structures/` are read-only copies of RCSB PDB mmCIF entries verified against `structures/MANIFEST.json` and `SHA256SUMS.txt`:

| PDB ID | Organism | Resolution / Method | Role in Workbench | Key Measured / Deposited Properties |
| :--- | :--- | :--- | :--- | :--- |
| **`6TAP`** | *Drosophila melanogaster* | 3.5 A (Cryo-EM) | Primary global capsid + radical-scavenging geometry | 60 icosahedral operators (`6TAP_operators.json`); 300 deposited polymer chains across two populations: 240 outer CA shell subunits (Chains `A-D`, $T=4$ surface lattice) + 60 interior NC zinc-knuckle fragments (Chain `E`, res 224-250) + 60 `Zn` ions (Chain `F`); $R_{\text{in}} = 15.54\text{ nm}$, $R_{\text{out}} = 19.93\text{ nm}$, $R_g = 17.55\text{ nm}$ |
| **`6TAR`** | *Drosophila melanogaster* | Cryo-EM | Local 5-fold capsomer reference | dArc1 five-fold capsomer |
| **`6S7X`** | *Drosophila melanogaster* | X-ray | Local CA dimer interface | dArc1 CA domain dimer (orthorhombic) |
| **`6S7Y`** | *Drosophila melanogaster* | X-ray | Local CA dimer interface | dArc1 CA domain dimer (hexagonal) |
| **`6GSE`** | *Rattus norvegicus* | Solution NMR ($C_1$) | Mammalian domain monomer | Rat Arc capsid-domain monomer (single chain; no icosahedral lattice) |
| **`4X3H`** | *Rattus norvegicus* | 2.40 A (X-ray) | Mammalian N-lobe domain | Rat Arc N-lobe bound to Stargazin peptide (domain only; no icosahedral lattice) |
| **`6R24`** | *Saccharomyces cerevisiae* | 7.5 A (Cryo-EM) | Retrotransposon global shell | Ty3/gypsy retrotransposon capsid (540 subunits, $T=9$, outer diameter $55.11\text{ nm}$, $R_g = 24.87\text{ nm}$) |
| **`1OHF`** | *Nudaurelia capensis* $\omega$ virus | X-ray | Canonical $T=4$ virus shell | NwV $T=4$ ssRNA capsid (240 subunits, outer diameter $41.79\text{ nm}$, $R_g = 16.89\text{ nm}$) |
| **`3H47`** | HIV-1 | X-ray | Retroviral CA hexamer | HIV-1 CA hexamer local interface reference |

> [!IMPORTANT]
> 1. **Domain-only vs Capsid Accessions**: Neither `6GSE` (rat Arc NMR monomer) nor `4X3H` (rat Arc N-lobe crystal domain) contains an icosahedral capsid lattice. All full-capsid icosahedral reconstructions use `6TAP` (*Drosophila* dArc1) alongside `6R24` (yeast Ty3/gypsy) and `1OHF` (NwV $T=4$).
> 2. **Two-Population Subunit Decomposition in `6TAP`**: $T=4$ ($60 \times 4 = 240$) applies strictly to the outer CA shell (`Chains A-D`). The `300` total deposited polymer chains sum two populations (`240` outer CA shell + `60` interior `Chain E` NC zinc-knuckle fragments, plus `60` `Chain F` $\text{Zn}^{2+}$ ions); `300/60 = 5` is not a Caspar-Klug triangulation number.
> 3. **Physical Scale & Regime (Radical Chemistry, Not Photon Shielding)**: Across a $39.85\text{ nm}$ protein shell at $\rho \approx 1.35\text{ g/cm}^3$ (mass thickness $\approx 5.4 \times 10^{-6}\text{ g/cm}^2$), $100\text{ keV}$ photon interaction probability is $\sim 10^{-6}$ (the shell is transparent to photons). The physically active regime at the $1\text{ to }50\text{ nm}$ scale is downstream of primary ionization: low-energy secondary electron track structure ($\sim 1\text{ keV}$ range $\approx 50\text{ nm}$ in water) and aqueous radical diffusion/scavenging ($\text{OH}^\bullet$ interception by aromatic/sulfur residues in the `240`-subunit CA shell and redox chemistry at the `60` interior $\text{Cys}_2\text{HisCys-Zn}^{2+}$ sites).

---

## 2. Four-Platform Visualization & Edge-Case Workflows

### A. UCSF ChimeraX 1.12 (`workbench/`)

- **Interactive Camera Rack on Verified Icosahedral Axes**:
  ```bash
  /Applications/ChimeraX-1.12.app/Contents/MacOS/ChimeraX workbench/arc_camera_rack.cxc
  ```
- **Portable 3-Scene AHRT Generator (regenerates `.cxs` sessions + `1920x1080` PNGs)**:
  ```bash
  /Applications/ChimeraX-1.12.app/Contents/MacOS/ChimeraX --nogui --offscreen workbench/scene_ahrt_01_portable.py
  ```
- **Pre-rendered `1920x1080` Outputs (`workbench/renders/`)**:
  - `workbench/renders/arc_structural_reference.png`: Mammalian `6GSE` monomer vs `4X3H` N-lobe domain.
  - `workbench/renders/arc_heterarchy.png`: Four-structure ASU comparison (`6GSE`, `6TAP`, `6R24`, `1OHF`).
  - `workbench/renders/arc_radiation_ready.png`: Full 60-operator `6TAP` capsid with front-hemisphere transparency revealing the 60 interior NC zinc-knuckle chains (`#ec4899`) and 60 `Zn` ions (`#facc15`).

### B. PyMOL (`pymol/`)

- **Headless Ray-Traced `1920x1080` Export**:
  ```bash
  /Applications/PyMOL.app/Contents/MacOS/PyMOL -cq pymol/arc_ahrt_capsid_scene.py
  ```
- **Interactive GUI Session**:
  ```bash
  /Applications/PyMOL.app/Contents/MacOS/PyMOL pymol/arc_ahrt_capsid_scene.pml
  ```
- **Pre-rendered `1920x1080` Outputs (`pymol/renders/`)**:
  - `pymol/renders/pymol_arc_domain_and_asu_comparison.png`: 4-panel grid view of `6GSE`, `4X3H`, `6TAP` ASU, and `6R24` ASU.
  - `pymol/renders/pymol_6tap_icosahedral_radiation_target.png`: Full 60-operator `6TAP` icosahedral assembly in PyMOL with cutaway front CA shell (`skyblue`, 82% transparent), back CA shell (`deepteal`), 60 interior NC zinc-knuckle chains (`hotpink`), and 60 `Zn` spheres (`tv_yellow`).
  - `pymol/renders/pymol_render_receipt.json`: SHA-256 and dimension receipt.

### C. MATLAB R2026b (`matlab/`)

- **4-Panel 3D Capsid + Level 3 Interface + TOPAS-nBio Radiolysis Dashboard**:
  ```bash
  /Applications/MATLAB_R2026b.app/Contents/MacOS/matlab -sd "$PWD" -batch "run(fullfile(pwd, 'matlab', 'arc_ahrt_3d_visualization.m'))"
  ```
- **V4-ASRA Scavenger Curve & Table Generator**:
  ```bash
  /Applications/MATLAB_R2026b.app/Contents/MacOS/matlab -sd "$PWD" -batch "run(fullfile(pwd, 'matlab', 'v4_asra_analysis.m'))"
  ```
- **Outputs**:
  - `matlab/renders/matlab_arc_ahrt_3d_dashboard.png` and `matlab/renders/matlab_arc_ahrt_3d_receipt.json`
  - `docs/v4_matlab_curves.png`, `docs/v4_matlab_series.csv`, and `docs/v4_matlab_receipt.json`

### D. RStudio / R (`rstudio/`)

- **Edge-Case Statistical & Physical Audit (`rstudio/arc_ahrt_edge_cases.R`)**:
  ```bash
  Rscript rstudio/arc_ahrt_edge_cases.R
  ```
  Or in RStudio:
  ```r
  source("rstudio/arc_ahrt_edge_cases.R")
  ```
- **Four Audited Edge Cases**:
  1. **EC1 (Symmetry-Only Interfaces & `A-B` Percolation Onset)**: In `geometry/interface_residue_contacts_sweep.json`, `A-A` ($C_5$ pentamer), `B-D`, and `C-C` ($C_3$ hexamer) have `0` contacts in the single ASU across all cutoffs (`3.5` to `6.0 A`) and only emerge under icosahedral operators (`83` combined unique pairs at `5.0 A`); `A-B` has `0` contacts at `<= 4.5 A` and emerges at `5.0 A` via a single inter-operator contact (`LYS181-ASN204` at `4.6059 A`, operator `2`).
  2. **EC2 (Level 3 Omnibus-12 Multiplicity Boundary)**: In `artifacts/m7_arc_level3_interface_comparison_receipt.json`, `8/12` permutation tests survive omnibus Benjamini-Hochberg FDR ($q \le 0.05$) while `3/12` survive conservative omnibus-12 Holm-Bonferroni ($p \le 0.05$); `homologous_83_residue_same_domain` passes within-family Holm ($p = 0.0325$) and omnibus BH ($q = 0.0195$) but sits on the omnibus-12 Holm boundary ($p = 0.0585$).
  3. **EC3 (V4-ASRA Radical Floor & Secondary Radical Sparing)**: In `docs/v4_matlab_series.csv`, `W-C1` (`25 mM NO3-` + `MeOH`) drives $G(\text{OH}^0)$ to an exact `0` floor at $\ge 0.1\text{ M}$ while $G(\text{H}^0)$ rises by $+67.9\%$ (`0.3510 -> 0.5893` molecules/100 eV); in `W-C2` (`O2`-saturated `DMSO`), $G(e_{\text{aq}}^{-1})$ rises by $+23.6\%$ (`0.01543 -> 0.01908`) as $\text{OH}^0$ scavenging suppresses $\text{OH} + e_{\text{aq}}^-$ recombination.
  4. **EC4 (Chain `E` NC Zinc-Knuckle B-Factor Dispersion)**: In `6TAP`, interior Chain `E` (`residues 224-250`) has a B-factor standard deviation of `32.25 A^2` (`range 9.96-160.89 A^2`), more than $3\times$ the dispersion of the outer CA shell chains `A-D` (`7.68-10.61 A^2`).
- **Outputs (`rstudio/renders/`)**:
  - `rstudio/renders/rstudio_arc_ahrt_edge_cases.png` (`1920x1080`)
  - `rstudio/renders/rstudio_edge_cases_summary.csv`
  - `rstudio/renders/rstudio_edge_cases_receipt.json`

---

## 3. Reproducibility & Verification Suite

Run the integrity and unit test checks from the repository root:

```bash
# 1. Verify all 22 SHA-256 pinned baseline structural and Level 1-3 receipts
shasum -a 256 -c SHA256SUMS.txt

# 2. Run M7 ARC Level 2 & Level 3 interface graph & permutation test suites (22 tests)
pytest code/test_m7_arc_level2_interface_atlas.py code/test_m7_arc_level3_interface_comparison.py -q

# 3. Run AHRT Stage A photon transport (NIST XCOM + Klein-Nishina) unit tests (11 tests)
PYTHONPATH=src python3 -m unittest discover -s tests -v
```
