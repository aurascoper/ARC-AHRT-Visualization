% ==============================================================================
% matlab/v4_asra_analysis.m
% AHRT V4-ASRA MATLAB Analysis of TOPAS-nBio W-C1 and W-C2 Scavenger Runs
% Portable paths via mfilename('fullpath'). Strictly zero em-dashes or en-dashes.
% ==============================================================================

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = fullfile(pwd, 'matlab');
end
repoRoot = fileparts(scriptDir);
rawPath = fullfile(repoRoot, 'data', 'wa04_raw.json');
docsDir = fullfile(repoRoot, 'docs');

raw = jsondecode(fileread(rawPath));
runs = raw.runs;
fn = fieldnames(runs);

species = {'H2O2^0', 'OH^0', 'e_aq^-1', 'H_2^0', 'H^0'};
rows = {};
for i = 1:numel(fn)
    k = fn{i};
    r = runs.(k);
    st = r.completion_state;
    tok = regexp(k, '_([0-9.eE+-]+)M$', 'tokens');
    if isempty(tok), continue; end
    conc = str2double(tok{1}{1});
    if startsWith(k, 'wc1_'), fix = 1; else, fix = 2; end
    for s = 1:numel(species)
        sp = species{s};
        sp_field = matlab.lang.makeValidName(sp);
        if isfield(r.final_G, sp_field)
            g = r.final_G.(sp_field).G;
            sd = r.final_G.(sp_field).err_1sd;
        else
            flds = fieldnames(r.final_G);
            g = NaN; sd = NaN;
            for ff = 1:numel(flds)
                if strcmp(regexprep(sp, '[^a-zA-Z0-9]', '_'), regexprep(flds{ff}, '[^a-zA-Z0-9]', '_'))
                    g = r.final_G.(flds{ff}).G;
                    sd = r.final_G.(flds{ff}).err_1sd;
                end
            end
        end
        rows(end+1, :) = {k, conc, fix, sp, g, sd, st}; %#ok<AGROW>
    end
end

T = cell2table(rows, 'VariableNames', {'config', 'conc', 'fixture', 'species', 'G_1us', 'err1sd', 'state'});
T = sortrows(T, {'fixture', 'conc', 'species'});
writetable(T, fullfile(docsDir, 'v4_matlab_series.csv'));

f = figure('Visible', 'off', 'Position', [100 100 1920 1080], 'Color', 'w');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for fix = 1:2
    nexttile; hold on; set(gca, 'XScale', 'log', 'FontSize', 12);
    for s = 1:numel(species)
        sp = species{s};
        sub = T(T.fixture == fix & strcmp(T.species, sp) & strcmp(T.state, 'COMPLETED'), :);
        if ~isempty(sub)
            errorbar(sub.conc, sub.G_1us, sub.err1sd, '-o', 'LineWidth', 1.8, 'DisplayName', sp);
        end
    end
    grid on; xlabel('Scavenger concentration (M)'); ylabel('G at 1 us (molecules / 100 eV)');
    if fix == 1
        title('W-C1: 25 mM NO3- + MeOH (COMPLETED runs)');
    else
        title('W-C2: O2-sat + DMSO (COMPLETED runs)');
    end
    legend('Location', 'best');
end

exportgraphics(f, fullfile(docsDir, 'v4_matlab_curves.png'), 'Resolution', 200);
close(f);

rec = struct('matlab_version', version, 'total_configs', numel(fn), ...
    'completed_configs', sum(strcmp(T.state, 'COMPLETED')) / numel(species), ...
    'csv_rows', height(T), 'timestamp_utc', char(datetime('now', 'TimeZone', 'UTC', 'Format', 'yyyy-MM-dd''T''HH:mm:ss''Z''')));
fid = fopen(fullfile(docsDir, 'v4_matlab_receipt.json'), 'w');
fprintf(fid, '%s\n', jsonencode(rec, 'PrettyPrint', true));
fclose(fid);
disp('MATLAB V4-ASRA analysis complete.');
