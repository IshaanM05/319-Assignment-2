function headtohead_figure(R, suite_names, suite_fields, figdir)
%HEADTOHEAD_FIGURE  Draw the ablation figures from the aggregated results in R.
%   R.(suite).{config} holds the mean metrics produced by run_headtohead.m.
%   Writes the compact two-panel figure used in the paper
%   (figures/fig_headtohead.{png,pdf}) and a full six-panel diagnostic version
%   (figures/fig_headtohead_full.{png,pdf}).  Kept in a separate function so the
%   figures can be redrawn from results/headtohead.mat without re-simulating.
sS = 'Stress24';  sO = 'Orig40';
getv = @(S,fn) cellfun(@(r) r.(fn), R.(S));
barlab = {'Base','+FF','Final','aQ','IMM'};
ns = numel(suite_names);

% ---- compact two-panel version (paper figure) ------------------------------
f = figure('Position',[50 50 760 300],'Visible','off');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
ms = zeros(numel(barlab), 2);  th = ms;
for q = 1:ns
    ms(:,q) = getv(suite_fields{q},'rms_miss');
    th(:,q) = getv(suite_fields{q},'first_hit');
end
nexttile; bar(ms); grid on; set(gca,'xticklabel',barlab,'FontSize',7);
ylabel('RMS miss [m]'); title('Engagement accuracy'); ylim([0 1.6]);
legend(suite_names,'Location','northwest','FontSize',7);
nexttile; bar(th); grid on; set(gca,'xticklabel',barlab,'FontSize',7);
ylabel('interception time [s]'); title('Interception time'); ylim([0 4.4]);
exportgraphics(f, fullfile(figdir,'fig_headtohead.png'), 'Resolution', 300, 'BackgroundColor', 'white');
exportgraphics(f, fullfile(figdir,'fig_headtohead.pdf'), 'BackgroundColor', 'white');
close(f);

% ---- full six-panel diagnostic version ------------------------------------
f = figure('Position',[50 50 1100 620],'Visible','off');
tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile; bar(getv(sS,'elos_hirate_deg')); set(gca,'xticklabel',barlab); grid on;
 ylabel('deg'); title('Pointing error at high rate (stress)');
nexttile; bar(getv(sS,'max_elos_deg')); set(gca,'xticklabel',barlab); grid on;
 ylabel('deg'); title('Peak pointing error (stress)');
nexttile; bar(getv(sS,'rms_pred')); set(gca,'xticklabel',barlab); grid on;
 ylabel('m'); title('RMS prediction error (stress)');
nexttile; bar(getv(sS,'first_hit')); set(gca,'xticklabel',barlab); grid on;
 ylabel('s'); title('Interception time (stress)');
nexttile;
 sr = [getv(sO,'success_rate'); getv(sS,'success_rate')]';
 bar(sr); set(gca,'xticklabel',barlab); grid on; ylabel('%'); ylim([0 105]);
 title('Interception success rate'); legend(suite_names,'Location','south','FontSize',7);
nexttile; bar(getv(sO,'elos_acq_deg')); set(gca,'xticklabel',barlab); grid on;
 ylabel('deg'); title('RMS pointing error (original 40)');
exportgraphics(f, fullfile(figdir,'fig_headtohead_full.png'), 'Resolution', 200, 'BackgroundColor', 'white');
exportgraphics(f, fullfile(figdir,'fig_headtohead_full.pdf'), 'BackgroundColor', 'white');
close(f);
end
