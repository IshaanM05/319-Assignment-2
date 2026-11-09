%RUN_HEADTOHEAD  Head-to-head benchmark: feedback-only baseline -> finalized design.
%
%   Compares, with ONE frozen design per column (no per-case retuning):
%     - Baseline (lag)    : lag-matched LTI feedback only, no feedforward, no dual mode
%     - + feedforward     : + the model-inverse 2-DOF feedforward (LTI)
%     - Finalized (FF+DM) : + the bounded dual-mode acquisition front-end  == DEFAULT
%     - final/est=aQ      : finalized, estimator swapped to NIS-adaptive Q
%     - final/est=IMM     : finalized, estimator swapped to CV/CA/CT IMM
%
%   over two suites:
%     - the ORIGINAL 40 randomized seeds (101-140): the finalized design must not
%       regress here
%     - the STRESS suite (sustained high-rate coordinated turns with a full 0.5 s
%       measurement loss inside the turn + a sharp tangential kick): the regime
%       that exposes the CA filter's and the lag match's admitted limits.
%
%   Writes results/headtohead.csv, results/headtohead.tex,
%   figures/fig_headtohead.{png,pdf}
clear; close all; clc;
paper_figure_style();
root = fileparts(mfilename('fullpath'));
figdir = fullfile(root,'figures'); resdir = fullfile(root,'results');
if ~isfolder(figdir), mkdir(figdir); end
if ~isfolder(resdir), mkdir(resdir); end

% Ablation of the FINALIZED design (= defaults: CA + lag-matched feedback +
% feedforward + dual-mode).  Each column is one frozen design, no per-case retuning.
cfgs  = {struct('ctrl','lag','acq','linear'), ...   % LTI lag feedback only (prior-art baseline)
         struct('acq','linear'), ...                % + feedforward (LTI 2-DOF)
         struct(), ...                              % + dual-mode acquisition  == FINALIZED
         struct('est','adaptQ'), ...                % finalized, estimator swapped to adaptive-Q
         struct('est','IMM')};                      % finalized, estimator swapped to IMM
names = {'Baseline (lag)','+ feedforward','Finalized (FF+DM)','final/est=aQ','final/est=IMM'};
Nc = numel(cfgs);

seeds  = 101:140;                 % original robustness suite
Kstr   = 1:24;                    % stress suite
suite_names  = {'Original-40','Stress-24'};
suite_fields = {'Orig40','Stress24'};
suite_cases  = { arrayfun(@(s) struct('seed',s),  seeds, 'uni',0), ...
                 arrayfun(@(k) struct('stress',k), Kstr,  'uni',0) };

% metric fields pulled from L.metrics (elevation errors converted to deg later)
flds = {'success','t_first_hit','rms_pred','rms_pred_loss','rms_pos','rms_pos_loss', ...
        'rms_elos_acq','rms_elos_hirate','max_elos_acq','rms_miss_acq_buf','t_acq','Vsat_frac'};

R = struct();   % R.(suitefield){cfg} = aggregated row
fprintf('Running head-to-head: %d configs x (%d + %d) cases ...\n', Nc, numel(seeds), numel(Kstr));
tic;
for si = 1:numel(suite_names)
    sname = suite_fields{si};
    cases = suite_cases{si};
    for c = 1:Nc
        acc = struct(); for f = 1:numel(flds), acc.(flds{f}) = nan(numel(cases),1); end
        for ic = 1:numel(cases)
            L = simulate_case(cases{ic}, cfgs{c});  m = L.metrics;
            for f = 1:numel(flds), acc.(flds{f})(ic) = m.(flds{f}); end
        end
        row.success_rate   = 100*mean(acc.success);
        row.first_hit      = mean(acc.t_first_hit,'omitnan');
        row.rms_pred       = mean(acc.rms_pred);
        row.rms_pred_loss  = mean(acc.rms_pred_loss,'omitnan');
        row.rms_pos        = mean(acc.rms_pos);
        row.rms_pos_loss   = mean(acc.rms_pos_loss,'omitnan');
        row.elos_acq_deg   = rad2deg(mean(acc.rms_elos_acq,'omitnan'));
        row.elos_hirate_deg= rad2deg(mean(acc.rms_elos_hirate,'omitnan'));
        row.max_elos_deg   = rad2deg(mean(acc.max_elos_acq,'omitnan'));
        row.rms_miss       = mean(acc.rms_miss_acq_buf,'omitnan');
        row.t_acq          = mean(acc.t_acq,'omitnan');
        row.vsat_pct       = 100*mean(acc.Vsat_frac);
        R.(sname){c} = row;
    end
end
fprintf('Done in %.1f s\n', toc);

%% ---- console + CSV ----------------------------------------------------------
mfields = {'success_rate','first_hit','rms_pred','rms_pred_loss','rms_pos','rms_pos_loss', ...
           'elos_acq_deg','elos_hirate_deg','max_elos_deg','rms_miss','t_acq','vsat_pct'};
hdr = [{'suite','config'}, mfields];
rows = {};
for si = 1:numel(suite_names)
    sname = suite_fields{si};
    fprintf('\n=== %s ===\n', suite_names{si});
    fprintf('%-16s %6s %7s %7s %8s %7s %8s %8s %9s %8s\n', 'config','succ%','tHit','pred','predLos','posLos','elosAcq','hiRate','maxElos','miss');
    for c = 1:Nc
        r = R.(sname){c};
        fprintf('%-16s %6.1f %7.2f %7.2f %8.2f %7.2f %8.4f %8.4f %9.4f %8.2f\n', names{c}, ...
            r.success_rate, r.first_hit, r.rms_pred, r.rms_pred_loss, r.rms_pos_loss, ...
            r.elos_acq_deg, r.elos_hirate_deg, r.max_elos_deg, r.rms_miss);
        rowc = {suite_names{si}, names{c}};
        for f = 1:numel(mfields), rowc{end+1} = r.(mfields{f}); end %#ok<SAGROW>
        rows(end+1,:) = rowc; %#ok<SAGROW>
    end
end
writecell([hdr; rows], fullfile(resdir,'headtohead.csv'));
save(fullfile(resdir,'headtohead.mat'), 'R', 'names', 'suite_names', 'suite_fields');

%% ---- LaTeX table (stress suite: where the upgrades earn their keep) -------------
sname = 'Stress24';
fid = fopen(fullfile(resdir,'headtohead.tex'),'w');
fprintf(fid, '%% Auto-generated by run_headtohead.m -- stress suite (24 cases), one frozen design/column\n');
fprintf(fid, '\\begin{tabular}{@{}lccccc@{}}\n\\toprule\n');
fprintf(fid, ' & Baseline & +FF & \\textbf{Final} & est=aQ & est=IMM \\\\\n\\midrule\n');
pr = @(fn,fmt) fprintf(fid, fmt, R.(sname){1}.(fn), R.(sname){2}.(fn), R.(sname){3}.(fn), R.(sname){4}.(fn), R.(sname){5}.(fn));
fprintf(fid,'Success rate [\\%%]'); pr('success_rate',' & %.1f & %.1f & %.1f & %.1f & %.1f \\\\\n');
fprintf(fid,'Interception time [s]'); pr('first_hit',' & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n');
fprintf(fid,'RMS pred.\\ error [m]'); pr('rms_pred',' & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n');
fprintf(fid,'RMS pred.\\ (coast) [m]'); pr('rms_pred_loss',' & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n');
fprintf(fid,'Pointing at high rate [deg]'); pr('elos_hirate_deg',' & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n');
fprintf(fid,'Peak pointing [deg]'); pr('max_elos_deg',' & %.3f & %.3f & %.3f & %.3f & %.3f \\\\\n');
fprintf(fid,'RMS miss [m]'); pr('rms_miss',' & %.2f & %.2f & %.2f & %.2f & %.2f \\\\\n');
fprintf(fid, '\\bottomrule\n\\end{tabular}\n');
fclose(fid);

%% ---- figures ---------------------------------------------------------------
headtohead_figure(R, suite_names, suite_fields, figdir);
disp('run_headtohead done.');