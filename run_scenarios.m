%RUN_SCENARIOS  Standardized scenarios A, B, C over the full 0..10 s: figures + metrics table.
%   Writes figures/fig_scenario_{A,B,C}.{png,pdf}, figures/fig_cov.{png,pdf},
%   results/scenario_metrics.csv, results/scenario_metrics.tex, results/scenarios.mat
clear; close all; clc;
paper_figure_style();
root = fileparts(mfilename('fullpath'));
figdir = fullfile(root,'figures'); resdir = fullfile(root,'results');
if ~isfolder(figdir), mkdir(figdir); end
if ~isfolder(resdir), mkdir(resdir); end

cases = {'A','B','C', struct('base','B','loss',[4.0 0.5])};      % 4th: extra loss during lock-on (illustration)
names = {'A','B','C','C2'};
R = cell(1,numel(cases));
for i = 1:numel(cases)
    R{i} = simulate_case(cases{i});
    plot_scenario(R{i}, names{i}, figdir);
end
plot_cov(R, names, figdir);

%% metrics table ---------------------------------------------------------------------
hdr = {'case','RMS pos err [m]','RMS vel err [m/s]','RMS pred err [m]','RMS e_ref full [deg]','RMS e_los full [deg]', ...
       'RMS e_ref acq [deg]','RMS e_los acq [deg]','t_acq [s]','RMS miss (buf) [m]','max miss (buf) [m]','pulses','hits','first hit [s]','max pulse miss [m]','V sat [%]'};
T = cell(numel(R), numel(hdr));
fprintf('\n%-4s %7s %7s %7s | %8s %8s | %8s %8s %6s %8s %8s | %3s %3s %6s %7s %6s\n', 'case','posRMS','velRMS','prdRMS','eref','elos','eref_a','elos_a','tacq','rmsMissB','maxMissB','np','nh','t1hit','maxPls','Vsat%');
for i = 1:numel(R)
    M = R{i}.metrics;
    row = {names{i}, M.rms_pos, M.rms_vel, M.rms_pred, rad2deg(M.rms_eref_full), rad2deg(M.rms_elos_full), ...
           rad2deg(M.rms_eref_acq), rad2deg(M.rms_elos_acq), M.t_acq, M.rms_miss_acq_buf, M.max_miss_acq_buf, M.n_pulses, M.n_success, M.t_first_hit, M.max_miss_pulse, 100*M.Vsat_frac};
    T(i,:) = row;
    fprintf('%-4s %7.2f %7.2f %7.2f | %8.3f %8.3f | %8.4f %8.4f %6.2f %8.2f %8.2f | %3d %3d %6.2f %7.2f %6.1f\n', row{:});
end
writecell([hdr; T], fullfile(resdir,'scenario_metrics.csv'));
fid = fopen(fullfile(resdir,'scenario_metrics.tex'),'w');
for i = 1:numel(R)
    M = R{i}.metrics;
    fprintf(fid, '%s & %.2f & %.2f & %.2f & %.3f & %.3f & %.2f & %d/%d & %.2f & %.2f \\\\\n', names{i}, M.rms_pos, M.rms_vel, M.rms_pred, ...
        rad2deg(M.rms_eref_acq), rad2deg(M.rms_elos_acq), M.t_acq, M.n_success, M.n_pulses, M.t_first_hit, M.max_miss_pulse);
end
fclose(fid);
save(fullfile(resdir,'scenarios.mat'), 'R', 'names');
% pulse tables
for i = 1:numel(R)
    fprintf('\nPulses, case %s  [t, true miss d_perp (m), success, range (m)]\n', names{i}); disp(R{i}.pulses);
end
% covariance numbers for the report (Scenario C)
Lc = R{3};  k0 = find(abs(Lc.t-1.0)<1e-9); k1 = find(abs(Lc.t-1.5)<1e-9);
fprintf('\nScenario C covariance: sigma_x,sigma_y at t=0.95 s : %.3f, %.3f m ; at t=1.49 s : %.3f, %.3f m ; at t=1.55 s : %.3f, %.3f m\n', ...
    Lc.sx(k0-50), Lc.sy(k0-50), Lc.sx(k1-10), Lc.sy(k1-10), Lc.sx(k1+50), Lc.sy(k1+50));
fprintf('                       sigma_vx,sigma_vy at t=0.95 : %.3f, %.3f ; 1.49 : %.3f, %.3f ; 1.55 : %.3f, %.3f m/s\n', ...
    Lc.svx(k0-50), Lc.svy(k0-50), Lc.svx(k1-10), Lc.svy(k1-10), Lc.svx(k1+50), Lc.svy(k1+50));
Lb = R{4};  k0 = find(abs(Lb.t-4.0)<1e-9); k1 = find(abs(Lb.t-4.5)<1e-9);
fprintf('Case C2 (loss 4.0-4.5 s): sigma_x,sigma_y at 3.95 : %.3f, %.3f ; 4.49 : %.3f, %.3f ; 4.55 : %.3f, %.3f m\n', ...
    Lb.sx(k0-50), Lb.sy(k0-50), Lb.sx(k1-10), Lb.sy(k1-10), Lb.sx(k1+50), Lb.sy(k1+50));
fprintf('  peak |theta_ref - theta_los| (deg) in C2 during loss+0.5 s: %.4f ; peak true miss d_perp in the same window: %.2f m\n', ...
    max(abs(rad2deg(Lb.thref(k0:k1+500) - Lb.theta_los(k0:k1+500)))), max(Lb.dperp(k0:k1+500)));

% =========================================================================
function shade_loss(ax, L)
yl = ylim(ax); d = diff([0; L.loss; 0]); a = find(d==1); b = find(d==-1);
for q = 1:numel(a)
    patch(ax, [L.t(a(q)) L.t(min(b(q),numel(L.t))) L.t(min(b(q),numel(L.t))) L.t(a(q))], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.9 0.2], 'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
ylim(ax, yl); uistack(findobj(ax,'Type','patch'),'bottom');
end

function plot_scenario(L, name, figdir)
f = figure('Position',[50 50 1000 800],'Visible','off','DefaultAxesFontSize',9);
tl = tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
d2 = @(v) rad2deg(v);
% (a) trajectory
ax = nexttile; plot(L.x, L.y, 'k', 'LineWidth', 1.6); hold on; plot(L.xh, L.yh, 'b--', 'LineWidth', 1.1); plot(L.xpx, L.xpy, 'r:', 'LineWidth', 1.3);
plot(L.x(1), L.y(1), 'ko', 'MarkerFaceColor', 'k'); axis equal; grid on; xlabel('x [m]'); ylabel('y [m]');
legend('true','estimated (now)','predicted (t+0.5 s)','Location','best','FontSize',7); title('(a) Target trajectory');
% (b) estimation / prediction errors
ax = nexttile; ex = L.xh - L.x; ey = L.yh - L.y;
np = round(L.P.tp/(L.t(2)-L.t(1))); m = 1:(numel(L.t)-np);
ep = hypot(L.xpx(m) - L.x(m+np), L.xpy(m) - L.y(m+np));
plot(L.t, ex, 'b', L.t, ey, 'Color', [0 0.6 0], 'LineWidth', 0.8); hold on; plot(L.t(m), ep, 'r', 'LineWidth', 1.0);
plot(L.t, 3*L.sx, 'b:', L.t, -3*L.sx, 'b:'); ylim([-10 15]); grid on; xlabel('t [s]'); ylabel('error [m]');
legend('x err','y err','|pred. err| (vs x(t+0.5))','\pm3\sigma_x','Location','northeast','FontSize',6); title('(b) Estimation and prediction error'); shade_loss(ax, L);
% (c) elevation
ax = nexttile; plot(L.t, d2(L.thref), 'r--', L.t, d2(L.th), 'b', L.t, d2(L.theta_los), 'k:', 'LineWidth', 1.3); grid on;
xlabel('t [s]'); ylabel('elevation [deg]'); legend('\theta_{ref}','\theta_{turret}','actual LOS','Location','southeast','FontSize',7); title('(c) Predicted reference, turret and LOS angle'); shade_loss(ax, L);
% (d) tracking error
ax = nexttile; plot(L.t, d2(L.thref - L.th), 'r', L.t, d2(L.theta_los - L.th), 'b', 'LineWidth', 1.0); grid on; ylim([-1.5 1.5]);
xlabel('t [s]'); ylabel('error [deg]'); legend('\theta_{ref}-\theta_{turret}','\theta_{LOS}-\theta_{turret}','Location','northeast','FontSize',7);
title('(d) Elevation errors (initial slew off-scale)'); shade_loss(ax, L);
% (e) miss distance and pulses
ax = nexttile; semilogy(L.t, max(L.dperp,1e-3), 'k', 'LineWidth', 0.8); hold on; yline(L.P.R_eng, 'r--'); grid on; ylim([0.01 3000]);
if ~isempty(L.pulses)
    ok = L.pulses(:,3) == 1;
    plot(L.pulses(ok,1), max(L.pulses(ok,2),1e-2), 'go', 'MarkerFaceColor', 'g'); plot(L.pulses(~ok,1), max(L.pulses(~ok,2),1e-2), 'rx', 'MarkerSize', 9, 'LineWidth', 2);
end
xlabel('t [s]'); ylabel('perp. miss distance [m]'); title('(e) Miss distance and engagement pulses (o hit, x miss), R_{eng}=2 m'); shade_loss(ax, L);
% (f) control voltage
ax = nexttile; plot(L.t, L.Vun, 'Color', [0.7 0.7 0.7]); hold on; plot(L.t, L.V, 'b', 'LineWidth', 1.2); yline(240, 'r:'); yline(-240, 'r:'); grid on; ylim([-300 300]);
xlabel('t [s]'); ylabel('V_m [V]'); legend('compensator state','V_m applied','Location','northeast','FontSize',7); title('(f) Motor voltage and saturation'); shade_loss(ax, L);
ttl = title(tl, sprintf('Scenario %s  (full 0-10 s;  yellow = measurement loss)', name));
ttl.Color = 'k';
exportgraphics(f, fullfile(figdir, ['fig_scenario_' name '.png']), 'Resolution', 250, 'BackgroundColor', 'white'); exportgraphics(f, fullfile(figdir, ['fig_scenario_' name '.pdf']), 'BackgroundColor', 'white');
close(f);
end

function plot_cov(R, names, figdir)
f = figure('Position',[50 50 900 320],'Visible','off');
cl = {'k','b','r','m'};
ax1 = subplot(1,2,1); hold on;
for i = [2 3 4], semilogy(R{i}.t, R{i}.sx, cl{i}, 'LineWidth', 1.1); end
set(gca,'YScale','log'); grid on; xlabel('t [s]'); ylabel('\sigma_x [m]'); xlim([0 6]);
legend(names{[2 3 4]}, 'Location', 'northeast', 'FontSize', 7); title('Position std. dev. \sigma_x (loss: C 1.0-1.5 s, C2 4.0-4.5 s)');
ax2 = subplot(1,2,2); hold on;
for i = [2 3 4], semilogy(R{i}.t, R{i}.sy, cl{i}, 'LineWidth', 1.1); end
set(gca,'YScale','log'); grid on; xlabel('t [s]'); ylabel('\sigma_y [m]'); xlim([0 6]); title('Position std. dev. \sigma_y');
exportgraphics(f, fullfile(figdir, 'fig_cov.png'), 'Resolution', 250, 'BackgroundColor', 'white'); exportgraphics(f, fullfile(figdir, 'fig_cov.pdf'), 'BackgroundColor', 'white'); close(f);
end
