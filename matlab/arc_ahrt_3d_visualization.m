% ==============================================================================
% matlab/arc_ahrt_3d_visualization.m
% Title: 4-Panel 1920x1080 3D Structural & Radiolysis Visualization Dashboard
% Epistemic Status: [MEASURED GEOMETRY & TOPAS-NBIO TRACK SERIES - NO MECHANISM CLAIMED]
% Compatible with MATLAB R2026b. Strictly zero em-dashes and zero en-dashes.
%
% Run from repository root:
%   /Applications/MATLAB_R2026b.app/Contents/MacOS/matlab -batch "run('matlab/arc_ahrt_3d_visualization.m')"
%
% Fail-closed: before drawing, the script checks the parsed structure against the
% Level 3 receipt's mmCIF inventory and checks the receipt against itself (the same
% identities scripts/render_6tap_interface_figure.py enforces). Any disagreement
% raises an error, so -batch exits nonzero and no figure or sidecar is written.
% Needs the JVM (the default) for SHA-256; do not run with -nojvm.
% ==============================================================================

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = fullfile(pwd, 'matlab');
end
repoRoot = fileparts(scriptDir);
outDir = fullfile(scriptDir, 'renders');
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% 1. Parse 6TAP.cif C-alpha (chains A-D, E) and ZN (chain F) + 60 operators
cifPath = fullfile(repoRoot, 'structures', '6TAP.cif');
opsPath = fullfile(repoRoot, 'structures', '6TAP_operators.json');
shellsPath = fullfile(repoRoot, 'geometry', 'shells.json');
l3Path = fullfile(repoRoot, 'artifacts', 'm7_arc_level3_interface_comparison_receipt.json');
v4CsvPath = fullfile(repoRoot, 'docs', 'v4_matlab_series.csv');

ops = jsondecode(fileread(opsPath));
shells = jsondecode(fileread(shellsPath));
l3 = jsondecode(fileread(l3Path));
v4Table = readtable(v4CsvPath);

% Extract ASU C-alpha and ZN coordinates from 6TAP.cif.
% Columns are located by their _atom_site header names, not by fixed positions.
cifLines = splitlines(fileread(cifPath));
hdrLines = strtrim(cifLines(startsWith(cifLines, '_atom_site.')));
col = @(name) find(strcmp(hdrLines, ['_atom_site.' name]), 1);
iAtom = col('label_atom_id'); iChain = col('label_asym_id');
iX = col('Cartn_x'); iY = col('Cartn_y'); iZ = col('Cartn_z');
if any(cellfun(@isempty, {iAtom, iChain, iX, iY, iZ}))
    error('arc:cifHeader', '6TAP.cif is missing one of the _atom_site columns this script reads');
end
atomMask = startsWith(cifLines, 'ATOM ') | startsWith(cifLines, 'HETATM ');
atomLines = cifLines(atomMask);

caShellASU = zeros(0, 3);
ncCoreASU = zeros(0, 3);
znASU = zeros(0, 3);

for i = 1:numel(atomLines)
    toks = strsplit(strtrim(atomLines{i}));
    if numel(toks) < numel(hdrLines)
        continue;
    end
    atomName = toks{iAtom};
    chainID = toks{iChain};
    xyz = str2double(toks([iX, iY, iZ]));
    if strcmp(atomName, 'CA') && ismember(chainID, {'A', 'B', 'C', 'D'})
        caShellASU(end+1, :) = xyz; %#ok<AGROW>
    elseif strcmp(atomName, 'CA') && strcmp(chainID, 'E')
        ncCoreASU(end+1, :) = xyz; %#ok<AGROW>
    elseif strcmp(chainID, 'F') || strcmp(atomName, 'ZN')
        znASU(end+1, :) = xyz; %#ok<AGROW>
    end
end

% ------------------------------------------------------------------------------
% Gates: parsed structure vs receipt inventory, and receipt vs itself
% ------------------------------------------------------------------------------
nOps = numel(ops);
prov = l3.mmcif_provenance;
if nOps ~= prov.inventory_reconciliation.operators_count
    error('arc:receiptMismatch', 'operators file has %d operators, receipt records %d', ...
        nOps, prov.inventory_reconciliation.operators_count);
end
inv = asCell(prov.struct_asym_inventory);
expected = struct('CA_CAPSID_SHELL', 0, 'NC_ZINC_KNUCKLE_FRAGMENT', 0, 'ZN_COORDINATION_LIGAND', 0);
for i = 1:numel(inv)
    role = inv{i}.functional_component_id;
    expected.(role) = expected.(role) + inv{i}.residue_count_per_asu * nOps;
end
parsed = [size(caShellASU, 1), size(ncCoreASU, 1), size(znASU, 1)] * nOps;
recorded = [expected.CA_CAPSID_SHELL, expected.NC_ZINC_KNUCKLE_FRAGMENT, expected.ZN_COORDINATION_LIGAND];
if ~isequal(parsed, recorded)
    error('arc:receiptMismatch', 'parsed CA/NC/Zn residues %s do not match receipt inventory %s', ...
        mat2str(parsed), mat2str(recorded));
end

profiles = asCell(l3.normalized_interface_chemistry_and_area_profiles.subinterface_profiles);
fingerprints = asCell(l3.symmetry_normalized_residue_fingerprints.subinterface_fingerprints);
fpNames = cellfun(@(f) f.interface_class, fingerprints, 'UniformOutput', false);
prNames = cellfun(@(p) p.interface_class, profiles, 'UniformOutput', false);
if ~isempty(setxor(fpNames, prNames))
    error('arc:receiptMismatch', 'receipt disagrees with itself: sub-interfaces in only one section: %s', ...
        strjoin(setxor(fpNames, prNames), ', '));
end
sumNN = 0; sumCC = 0; sumCN = 0;
for i = 1:numel(profiles)
    dap = profiles{i}.domain_architecture_profile;
    total = dap.NTD_NTD_pairs + dap.CTD_CTD_pairs + dap.NTD_CTD_pairs;
    fp = fingerprints{strcmp(fpNames, profiles{i}.interface_class)};
    if total ~= fp.unique_residue_pair_count_at_5A
        error('arc:receiptMismatch', 'receipt disagrees with itself: %s profile total %d ~= fingerprint pair count %d', ...
            profiles{i}.interface_class, total, fp.unique_residue_pair_count_at_5A);
    end
    sumNN = sumNN + dap.NTD_NTD_pairs; sumCC = sumCC + dap.CTD_CTD_pairs; sumCN = sumCN + dap.NTD_CTD_pairs;
end
audit = l3.constrained_graph_state_space_connectivity_and_mcmc_audit;
stub = audit.stub_difference_conservation_identity.empirical_6tap_5A_cutoff_values;
derived = [sumNN, sumCC, sumCN, sumNN + sumCC, sumCC - sumNN];
stubRec = [stub.observed_E_NN, stub.observed_E_CC, stub.observed_E_CN, stub.observed_E_same, stub.conserved_E_CC_minus_E_NN];
if ~isequal(derived, stubRec)
    error('arc:receiptMismatch', 'receipt disagrees with itself: summed rows give %s, stub-identity section records %s', ...
        mat2str(derived), mat2str(stubRec));
end
schemeObs = [audit.global_degree_preserving_simple_graph_double_edge_swap_mcmc.observed_same_domain_contacts, ...
    audit.within_subinterface_bipartite_degree_preserving_double_edge_swap_mcmc.observed_same_domain_contacts];
if any(schemeObs ~= sumNN + sumCC)
    error('arc:receiptMismatch', 'receipt disagrees with itself: null schemes record E_same %s, rows give %d', ...
        mat2str(schemeObs), sumNN + sumCC);
end
rOutNm = shells.x6TAP.descriptors.radial_profile.outer_radius_nm;  % p99 CA radius, from shells.json

% Subsample CA shell ASU by factor of 2 for crisp 60-copy 3D point cloud rendering
caShellSub = caShellASU(1:2:end, :);
allCA = zeros(size(caShellSub, 1) * nOps, 3);
allNC = zeros(size(ncCoreASU, 1) * nOps, 3);
allZN = zeros(size(znASU, 1) * nOps, 3);

nCA = size(caShellSub, 1);
nNC = size(ncCoreASU, 1);
nZN = size(znASU, 1);

for k = 1:nOps
    R = reshape(ops(k).matrix, [3, 3])';
    t = ops(k).vector(:)';
    allCA((k-1)*nCA + (1:nCA), :) = (caShellSub * R') + t;
    allNC((k-1)*nNC + (1:nNC), :) = (ncCoreASU * R') + t;
    allZN((k-1)*nZN + (1:nZN), :) = (znASU * R') + t;
end

% Convert Angstrom to nm and center at origin for 3D visualization
centerA = mean(allZN, 1);
allCA_nm = (allCA - centerA) / 10.0;
allNC_nm = (allNC - centerA) / 10.0;
allZN_nm = (allZN - centerA) / 10.0;

% Create 1920x1080 dark-theme figure
fig = figure('Visible', 'off', 'Position', [50, 50, 1920, 1080], 'Color', [0.03 0.05 0.08]);
tlo = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tlo, 'ARC Capsid Geometry & AHRT v0.1 Radiolysis Dashboard [MEASURED GEOMETRY - NO MECHANISM CLAIMED]', ...
    'Color', [0.97 0.98 0.99], 'FontSize', 15, 'FontWeight', 'bold');

% ------------------------------------------------------------------------------
% PANEL 1: 3D Cutaway of 60-Operator 6TAP Capsid (240 CA + 60 NC + 60 Zn)
% ------------------------------------------------------------------------------
ax1 = nexttile(tlo, 1);
hold(ax1, 'on');
set(ax1, 'Color', [0.05 0.08 0.12], 'XColor', [0.75 0.80 0.88], ...
    'YColor', [0.75 0.80 0.88], 'ZColor', [0.75 0.80 0.88], 'FontSize', 10);

% Cutaway front octant (x > 0 & y < 0 & z > 0, the one facing the camera) of CA shell so interior NC/Zn shows
cutMask = ~(allCA_nm(:,1) > 0 & allCA_nm(:,2) < 0 & allCA_nm(:,3) > 0);
scatter3(ax1, allCA_nm(cutMask,1), allCA_nm(cutMask,2), allCA_nm(cutMask,3), ...
    5, [0.22 0.74 0.97], 'filled', 'MarkerFaceAlpha', 0.35, 'DisplayName', 'CA Shell (Chains A-D, 240 copies)');
scatter3(ax1, allNC_nm(:,1), allNC_nm(:,2), allNC_nm(:,3), ...
    16, [0.93 0.28 0.60], 'filled', 'MarkerFaceAlpha', 0.85, 'DisplayName', 'NC Zn-Knuckle (Chain E, 60 copies)');
scatter3(ax1, allZN_nm(:,1), allZN_nm(:,2), allZN_nm(:,3), ...
    52, [0.98 0.80 0.08], 'filled', 'MarkerEdgeColor', [0.1 0.1 0.1], 'DisplayName', 'Zn2+ Ions (Chain F, 60 copies)');

axis(ax1, 'equal');
grid(ax1, 'on');
view(ax1, 38, 22);
xlabel(ax1, 'x (nm)'); ylabel(ax1, 'y (nm)'); zlabel(ax1, 'z (nm)');
title(ax1, sprintf('1. PDB 6TAP Full %d-Operator Assembly (Octant Cutaway, R_{out} = %.2f nm)', nOps, rOutNm), ...
    'Color', [0.97 0.98 0.99], 'FontSize', 12);
lg1 = legend(ax1, 'Location', 'southwest', 'TextColor', [0.95 0.96 0.98], 'Color', [0.08 0.11 0.16]);
lg1.FontSize = 9;

% ------------------------------------------------------------------------------
% PANEL 2: Topological Percolation Sweep (beta_0 & Algebraic Connectivity lambda_2)
% ------------------------------------------------------------------------------
ax2 = nexttile(tlo, 2);
hold(ax2, 'on');
set(ax2, 'Color', [0.05 0.08 0.12], 'XColor', [0.75 0.80 0.88], ...
    'YColor', [0.75 0.80 0.88], 'FontSize', 10);

sweep6TAP = shells.x6TAP.descriptors.contact_graph_sweep;
sweep6R24 = shells.x6R24.descriptors.contact_graph_sweep;
sweep1OHF = shells.x1OHF.descriptors.contact_graph_sweep;

cut_nm_6tap = arrayfun(@(s) s.cutoff_nm, sweep6TAP);
lam_6tap = arrayfun(@(s) s.algebraic_connectivity_dimensionless, sweep6TAP);
b0_6tap = arrayfun(@(s) s.beta_0_count, sweep6TAP);

cut_nm_6r24 = arrayfun(@(s) s.cutoff_nm, sweep6R24);
lam_6r24 = arrayfun(@(s) s.algebraic_connectivity_dimensionless, sweep6R24);

cut_nm_1ohf = arrayfun(@(s) s.cutoff_nm, sweep1OHF);
lam_1ohf = arrayfun(@(s) s.algebraic_connectivity_dimensionless, sweep1OHF);

yyaxis(ax2, 'left');
plot(ax2, cut_nm_6tap, lam_6tap, '-o', 'Color', [0.22 0.74 0.97], 'LineWidth', 2.0, 'DisplayName', '6TAP dArc1 \lambda_2');
plot(ax2, cut_nm_6r24, lam_6r24, '-s', 'Color', [0.98 0.45 0.09], 'LineWidth', 2.0, 'DisplayName', '6R24 Ty3/gypsy \lambda_2');
plot(ax2, cut_nm_1ohf, lam_1ohf, '-^', 'Color', [0.13 0.77 0.37], 'LineWidth', 2.0, 'DisplayName', '1OHF NwV T=4 \lambda_2');
ylabel(ax2, 'Algebraic Connectivity \lambda_2 (dimensionless)');
ax2.YAxis(1).Color = [0.22 0.74 0.97];

yyaxis(ax2, 'right');
plot(ax2, cut_nm_6tap, b0_6tap, '--d', 'Color', [0.98 0.80 0.08], 'LineWidth', 1.8, 'DisplayName', '6TAP Connected Components \beta_0');
ylabel(ax2, 'Connected Components \beta_0 (Int64 count)');
ax2.YAxis(2).Color = [0.98 0.80 0.08];

grid(ax2, 'on');
xlabel(ax2, 'Centroid Contact Cutoff (nm)');
title(ax2, '2. Capsid Contact-Graph Percolation Sweep (\lambda_2 & \beta_0 vs Cutoff)', ...
    'Color', [0.97 0.98 0.99], 'FontSize', 12);
lg2 = legend(ax2, 'Location', 'northwest', 'TextColor', [0.95 0.96 0.98], 'Color', [0.08 0.11 0.16]);
lg2.FontSize = 9;

% ------------------------------------------------------------------------------
% PANEL 3: Level 3 Subinterface BSA (Shrake-Rupley One-Side vs Voronoi Area)
% ------------------------------------------------------------------------------
ax3 = nexttile(tlo, 3);
hold(ax3, 'on');
set(ax3, 'Color', [0.05 0.08 0.12], 'XColor', [0.75 0.80 0.88], ...
    'YColor', [0.75 0.80 0.88], 'FontSize', 10);

nProf = numel(profiles);
labels = cell(nProf, 1);
srOneSide = zeros(nProf, 1);
vorArea = zeros(nProf, 1);

for i = 1:nProf
    p = profiles{i};
    labels{i} = strrep(p.interface_class, '_', ' ');
    srOneSide(i) = p.area_metrics_A2.shrake_rupley_per_side_bsa_A2;
    vorArea(i) = p.area_metrics_A2.voronoi_radical_plane_contact_area_A2;
end

b = bar(ax3, 1:nProf, [srOneSide, vorArea], 'grouped');
b(1).FaceColor = [0.22 0.74 0.97];
b(1).DisplayName = 'Shrake-Rupley Per-Side BSA (A^2)';
b(2).FaceColor = [0.66 0.33 0.97];
b(2).DisplayName = 'Voronoi Radical-Plane Area (A^2)';

set(ax3, 'XTick', 1:nProf, 'XTickLabel', labels, 'XTickLabelRotation', 28);
ylabel(ax3, 'Interface Area (A^2)');
grid(ax3, 'on');
title(ax3, '3. M7 ARC Level 3 Subinterface Area Calibration (10 Symmetry Classes)', ...
    'Color', [0.97 0.98 0.99], 'FontSize', 12);
lg3 = legend(ax3, 'Location', 'northeast', 'TextColor', [0.95 0.96 0.98], 'Color', [0.08 0.11 0.16]);
lg3.FontSize = 9;

% ------------------------------------------------------------------------------
% PANEL 4: AHRT V4-ASRA 1 us Radiolysis Scavenger Yields (W-C1 & W-C2)
% ------------------------------------------------------------------------------
ax4 = nexttile(tlo, 4);
hold(ax4, 'on');
set(ax4, 'Color', [0.05 0.08 0.12], 'XColor', [0.75 0.80 0.88], ...
    'YColor', [0.75 0.80 0.88], 'XScale', 'log', 'FontSize', 10);

wc1_oh = v4Table(v4Table.fixture == 1 & strcmp(v4Table.species, 'OH^0'), :);
wc1_h2o2 = v4Table(v4Table.fixture == 1 & strcmp(v4Table.species, 'H2O2^0'), :);
wc1_h = v4Table(v4Table.fixture == 1 & strcmp(v4Table.species, 'H^0'), :);
wc2_oh = v4Table(v4Table.fixture == 2 & strcmp(v4Table.species, 'OH^0'), :);
wc2_h2o2 = v4Table(v4Table.fixture == 2 & strcmp(v4Table.species, 'H2O2^0'), :);

errorbar(ax4, wc1_oh.conc, wc1_oh.G_1us, wc1_oh.err1sd, '-o', 'Color', [0.22 0.74 0.97], 'LineWidth', 1.8, 'DisplayName', 'W-C1 (MeOH): OH^0');
errorbar(ax4, wc1_h2o2.conc, wc1_h2o2.G_1us, wc1_h2o2.err1sd, '--s', 'Color', [0.13 0.77 0.37], 'LineWidth', 1.8, 'DisplayName', 'W-C1 (MeOH): H2O2^0');
errorbar(ax4, wc1_h.conc, wc1_h.G_1us, wc1_h.err1sd, ':d', 'Color', [0.98 0.80 0.08], 'LineWidth', 1.8, 'DisplayName', 'W-C1 (MeOH): H^0 (secondary rise)');
errorbar(ax4, wc2_oh.conc, wc2_oh.G_1us, wc2_oh.err1sd, '-^', 'Color', [0.98 0.45 0.09], 'LineWidth', 1.8, 'DisplayName', 'W-C2 (DMSO): OH^0');
errorbar(ax4, wc2_h2o2.conc, wc2_h2o2.G_1us, wc2_h2o2.err1sd, '--v', 'Color', [0.93 0.28 0.60], 'LineWidth', 1.8, 'DisplayName', 'W-C2 (DMSO): H2O2^0');

grid(ax4, 'on');
xlabel(ax4, 'Scavenger Concentration (M)');
ylabel(ax4, 'G at 1 us (molecules / 100 eV)');
title(ax4, '4. AHRT V4-ASRA TOPAS-nBio 1 us Scavenger Yields (W-C1 MeOH vs W-C2 DMSO)', ...
    'Color', [0.97 0.98 0.99], 'FontSize', 12);
lg4 = legend(ax4, 'Location', 'northeast', 'TextColor', [0.95 0.96 0.98], 'Color', [0.08 0.11 0.16]);
lg4.FontSize = 9;

outPng = fullfile(outDir, 'matlab_arc_ahrt_3d_dashboard.png');
exportgraphics(fig, outPng, 'Resolution', 150, 'BackgroundColor', [0.03 0.05 0.08]);
close(fig);

imgInfo = imfinfo(outPng);
inputs = struct( ...
    'level3_receipt', sha256File(l3Path), ...
    'mmcif_6tap', sha256File(cifPath), ...
    'operators', sha256File(opsPath), ...
    'shells', sha256File(shellsPath), ...
    'v4_series_csv', sha256File(v4CsvPath));
receipt = struct( ...
    'schema', 'arc-ahrt-matlab-3d-visualization-receipt-v2', ...
    'matlab_version', version, ...
    'figure', struct('file', 'matlab/renders/matlab_arc_ahrt_3d_dashboard.png', 'sha256', sha256File(outPng), ...
        'width_px', int64(imgInfo.Width), 'height_px', int64(imgInfo.Height)), ...
    'inputs_sha256', inputs, ...
    'reproducibility_note', ['Inputs are identified by sha256. The figure''s own sha256 reproduces only under ' ...
        'the same matlab_version; across versions, compare the numbers in this file instead.'], ...
    'operator_count', int64(nOps), ...
    'residues_in_assembly', struct('ca', int64(parsed(1)), 'nc', int64(parsed(2)), 'zn', int64(parsed(3))), ...
    'ca_points_rendered', int64(size(allCA, 1)), ...
    'ca_subsample_step', int64(2), ...
    'nc_points_rendered', int64(size(allNC, 1)), ...
    'zn_points_rendered', int64(size(allZN, 1)), ...
    'outer_radius_nm_from_shells_json', rOutNm, ...
    'receipt_identities', struct('E_NN', int64(sumNN), 'E_CC', int64(sumCC), 'E_CN', int64(sumCN), ...
        'E_same', int64(sumNN + sumCC), 'E_CC_minus_E_NN', int64(sumCC - sumNN)) ...
);
fid = fopen(fullfile(outDir, 'matlab_arc_ahrt_3d_receipt.json'), 'w');
fprintf(fid, '%s\n', jsonencode(receipt, 'PrettyPrint', true));
fclose(fid);
disp('MATLAB 3D ARC/AHRT visualization complete.');

% ------------------------------------------------------------------------------
% Local functions
% ------------------------------------------------------------------------------
function c = asCell(x)
% jsondecode returns a struct array for uniform JSON arrays and a cell array otherwise.
if iscell(x)
    c = x(:)';
else
    c = num2cell(x(:)');
end
end

function hex = sha256File(path)
fid = fopen(path, 'r');
if fid < 0
    error('arc:io', 'cannot open %s', path);
end
bytes = fread(fid, Inf, '*uint8');
fclose(fid);
md = java.security.MessageDigest.getInstance('SHA-256');
md.update(typecast(bytes, 'int8'));
hex = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2)', 1, []));
end
