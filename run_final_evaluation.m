%RUN_FINAL_EVALUATION  Robustness sweep across randomized bounded manoeuvres and
%   random loss intervals (the instructor has not released evaluation cases at the
%   time of writing, so Sec. 5/9 of the brief is satisfied with N_CASES >= 20
%   randomized trials, each with its own seed -> its own manoeuvre AND its own
%   loss interval(s), all within the stated envelope).  NO retuning: the exact
%   turret_controller.m / turret_params.m used for Scenarios A/B/C is reused.
%
%   Writes results/final_eval.csv, results/final_eval_summary.txt,
%   figures/fig_final_eval.{png,pdf}, results/final_eval.mat
clear; close all; clc;
root = fileparts(mfilename('fullpath'));
figdir = fullfile(root,'figures'); resdir = fullfile(root,'results');
if ~isfolder(figdir), mkdir(figdir); end
if ~isfolder(resdir), mkdir(resdir); end

N = 40;                                    % >= 20 required; 40 used for a tighter CI
seeds = 101:(100+N);
M = struct('succ',zeros(N,1),'thit',nan(N,1),'npulse',zeros(N,1),'nhit',zeros(N,1), ...
           'maxmiss',nan(N,1),'rmsmiss',nan(N,1),'firstmiss',nan(N,1),'rmspos',zeros(N,1),'rmsvel',zeros(N,1), ...
           'rmspred',zeros(N,1),'rmselos',zeros(N,1),'vsat',zeros(N,1),'tacq',nan(N,1),'nloss',zeros(N,1));

fprintf('Running %d randomized cases ...\n', N);
tAll = tic;
for i = 1:N
    L = simulate_case(struct('seed', seeds(i)));
    m = L.metrics;
    M.succ(i) = m.success; M.thit(i) = m.t_first_hit; M.npulse(i) = m.n_pulses; M.nhit(i) = m.n_success;
    M.maxmiss(i) = m.max_miss_acq_buf; M.rmsmiss(i) = m.rms_miss_acq_buf; M.firstmiss(i) = m.first_pulse_miss;
    M.rmspos(i) = m.rms_pos; M.rmsvel(i) = m.rms_vel; M.rmspred(i) = m.rms_pred;
    M.rmselos(i) = rad2deg(m.rms_elos_acq); M.vsat(i) = 100*m.Vsat_frac; M.tacq(i) = m.t_acq;
    M.nloss(i) = size(L.y(L.loss>0),1) > 0; %#ok<NASGU>
    M.nloss(i) = any(L.loss > 0);
    if mod(i,10)==0, fprintf('  %d/%d done (%.0f s)\n', i, N, toc(tAll)); end
end
fprintf('Total time: %.1f s\n', toc(tAll));

succ_rate = mean(M.succ);
fprintf('\n=== FINAL EVALUATION: %d randomized cases, no retuning ===\n', N);
fprintf('Interception success rate : %d/%d = %.1f %%\n', sum(M.succ), N, 100*succ_rate);
fprintf('First-hit time  (hits)    : mean %.3f s, median %.3f s, std %.3f s, min %.3f, max %.3f\n', ...
    mean(M.thit,'omitnan'), median(M.thit,'omitnan'), std(M.thit,'omitnan'), min(M.thit), max(M.thit));
fprintf('Pulses issued per case    : mean %.2f (range %d-%d)\n', mean(M.npulse), min(M.npulse), max(M.npulse));
fprintf('Successful-pulse fraction : %.1f %% of all issued pulses\n', 100*sum(M.nhit)/sum(M.npulse));
fprintf('Max miss distance (post-acquisition, true): mean %.2f m, 95th pct %.2f m, max %.2f m\n', ...
    mean(M.maxmiss), prctile(M.maxmiss,95), max(M.maxmiss));
fprintf('First-pulse true miss distance: mean %.2f m, median %.2f m, max %.2f m (<=2 m => hit)\n', ...
    mean(M.firstmiss), median(M.firstmiss), max(M.firstmiss));
fprintf('RMS position estimation error : mean %.2f m (range %.2f-%.2f)\n', mean(M.rmspos), min(M.rmspos), max(M.rmspos));
fprintf('RMS prediction error (0.5 s)   : mean %.2f m (range %.2f-%.2f)\n', mean(M.rmspred), min(M.rmspred), max(M.rmspred));
fprintf('RMS elevation error (post-acq) : mean %.4f deg (range %.4f-%.4f)\n', mean(M.rmselos), min(M.rmselos), max(M.rmselos));
fprintf('Acquisition time t_acq         : mean %.2f s, median %.2f s, max %.2f s\n', mean(M.tacq,'omitnan'), median(M.tacq,'omitnan'), max(M.tacq));
fprintf('Voltage-saturation duty cycle   : mean %.1f %% (range %.1f-%.1f)\n', mean(M.vsat), min(M.vsat), max(M.vsat));
fprintf('Cases that included a sensor-loss interval: %d/%d\n', sum(M.nloss), N);
if any(~M.succ)
    fprintf('Failed cases (seed): %s\n', mat2str(seeds(~M.succ)));
end

%% write CSV + summary ----------------------------------------------------------------
hdr = {'seed','success','t_first_hit','n_pulses','n_hits','max_miss_acq_buf_m','rms_miss_acq_buf_m','first_pulse_miss_m', ...
       'rms_pos_m','rms_vel_mps','rms_pred_m','rms_elos_acq_deg','t_acq_s','vsat_pct','had_loss'};
Tc = [seeds(:), M.succ, M.thit, M.npulse, M.nhit, M.maxmiss, M.rmsmiss, M.firstmiss, M.rmspos, M.rmsvel, M.rmspred, M.rmselos, M.tacq, M.vsat, M.nloss];
writecell([hdr; num2cell(Tc)], fullfile(resdir,'final_eval.csv'));

fid = fopen(fullfile(resdir,'final_eval_summary.txt'),'w');
fprintf(fid, 'FINAL EVALUATION: %d randomized cases (seeds %d-%d), ONE fixed design, no retuning\n', N, seeds(1), seeds(end));
fprintf(fid, 'Interception success rate: %d/%d = %.1f %%\n', sum(M.succ), N, 100*succ_rate);
fprintf(fid, 'First-hit time: mean %.3f s, median %.3f s, std %.3f s\n', mean(M.thit,'omitnan'), median(M.thit,'omitnan'), std(M.thit,'omitnan'));
fprintf(fid, 'Pulses/case: mean %.2f; successful-pulse fraction %.1f %%\n', mean(M.npulse), 100*sum(M.nhit)/sum(M.npulse));
fprintf(fid, 'Max miss distance (post-acq): mean %.2f m, 95th pct %.2f m, max %.2f m\n', mean(M.maxmiss), prctile(M.maxmiss,95), max(M.maxmiss));
fprintf(fid, 'RMS pos err: mean %.2f m; RMS pred err: mean %.2f m; RMS elos (post-acq): mean %.4f deg\n', mean(M.rmspos), mean(M.rmspred), mean(M.rmselos));
fprintf(fid, 'Acquisition time: mean %.2f s, median %.2f s, max %.2f s\n', mean(M.tacq,'omitnan'), median(M.tacq,'omitnan'), max(M.tacq));
fprintf(fid, 'Voltage saturation duty cycle: mean %.1f %%\n', mean(M.vsat));
fclose(fid);
save(fullfile(resdir,'final_eval.mat'), 'M', 'seeds', 'N');

%% figure: distributions ---------------------------------------------------------------
f = figure('Position',[50 50 1000 650],'Visible','off'); tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
nexttile; histogram(M.thit, 'BinWidth', 0.1, 'FaceColor', [0.2 0.4 0.8]); grid on; xlabel('first-hit time [s]'); ylabel('count');
title(sprintf('First-hit time (succ. %d/%d, mean %.2f s)', sum(M.succ), N, mean(M.thit,'omitnan')));
nexttile; histogram(M.rmsmiss, 'BinWidth', 0.2, 'FaceColor', [0.2 0.7 0.3]); grid on; xlabel('RMS miss dist., post-acq buffered [m]'); ylabel('count');
title('Steady-state tracking accuracy');
nexttile; histogram(M.rmspos, 'BinWidth', 0.1, 'FaceColor', [0.8 0.5 0.1]); grid on; xlabel('RMS position error [m]'); ylabel('count'); title('Estimator RMS error');
nexttile; histogram(M.rmspred, 'BinWidth', 0.2, 'FaceColor', [0.6 0.2 0.6]); grid on; xlabel('RMS prediction error [m]'); ylabel('count'); title('Predictor RMS error (0.5 s horizon)');
nexttile; histogram(M.npulse, 'BinMethod','integers', 'FaceColor', [0.3 0.3 0.3]); grid on; xlabel('pulses issued per case'); ylabel('count'); title('Pulses per 10 s engagement');
nexttile; scatter(seeds, M.thit, 28, M.succ, 'filled'); colormap(gca, [0.8 0.2 0.2; 0.2 0.7 0.2]); caxis([0 1]); grid on;
xlabel('case seed'); ylabel('first-hit time [s]'); title('Per-case outcome (green = hit, red = miss)');
exportgraphics(f, fullfile(figdir,'fig_final_eval.png'), 'Resolution', 250); exportgraphics(f, fullfile(figdir,'fig_final_eval.pdf'));
disp('run_final_evaluation done.');
