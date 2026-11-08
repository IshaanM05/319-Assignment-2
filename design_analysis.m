%DESIGN_ANALYSIS  Plant derivation check, compensator design, Bode/margins, closed-loop
%   analysis and anti-windup demonstration.  Reproduces Fig. 2-4 of the report and
%   writes results/design_numbers.txt.     Requires Control System Toolbox.
clear; close all; clc;
clear turret_params; P = turret_params();
root = fileparts(mfilename('fullpath'));
figdir = fullfile(root,'figures'); resdir = fullfile(root,'results');
if ~isfolder(figdir), mkdir(figdir); end
if ~isfolder(resdir), mkdir(resdir); end
fid = fopen(fullfile(resdir,'design_numbers.txt'),'w');
out = @(fmt, varargin) fprintf_both(fid, fmt, varargin{:});

%% 1. Plant from first principles ------------------------------------------------
s = tf('s');
Jeq = P.J_turret + P.N^2*P.J_m;
Gp  = P.N*P.Kt / ((Jeq*s^2 + P.b*s)*(P.Lm*s + P.Rm) + P.N^2*P.Kt*P.Ke*s);
% cross-check against the state-space model used in turret_plant.m
A = [0 1 0; 0 -P.b/Jeq P.N*P.Kt/Jeq; 0 -P.Ke*P.N/P.Lm -P.Rm/P.Lm];
B = [0; 0; 1/P.Lm];  Cm = [1 0 0];
Gss = tf(ss(A,B,Cm,0));
out('J_eq = %.1f kg m^2\n', Jeq);
out('Gp(s) = %s\n', evalc('zpk(Gp)'));
[num,den] = tfdata(Gp,'v'); k0 = num(end)/den(1);
out('Gp(s) = %.4f / [ s ( s^2 + %.3f s + %.2f ) ]\n', num(end)/den(1), den(2)/den(1), den(3)/den(1));
out('max |Gp - Gss| over 0.01..1e3 rad/s = %.2e\n', max(abs(squeeze(freqresp(Gp - Gss, logspace(-2,3,400))))));
pl = pole(Gp);  wn = abs(pl(2));  zeta = -real(pl(2))/wn;
out('poles: %s\n', mat2str(pl, 5));
out('resonance wn = %.2f rad/s, zeta = %.4f\n', wn, zeta);
Kdc = P.Kdc;  out('low-frequency velocity gain Kdc = N Kt/(N^2 Kt Ke + Rm b) = %.4e rad/s/V (=1/%.1f V s/rad)\n', Kdc, 1/Kdc);
out('max slew rate at Vmax = %.1f V: %.4f rad/s = %.2f deg/s\n', P.Vmax, P.Vmax*Kdc, rad2deg(P.Vmax*Kdc));
[~, ~] = deal(0);
% peak open-loop plant gain at resonance
out('|Gp(j wn)| = %.3e  (%.1f dB)\n', abs(freqresp(Gp, wn)), 20*log10(abs(freqresp(Gp, wn))));

%% 2. Compensators ----------------------------------------------------------------
Gc  = P.Kp_eff/(1 + s/P.p_eff);                    % FINAL
Gpi = P.Kp_pi*(1 + 1/(P.Ti_pi*s));                 % suggested PI (for comparison)
L   = Gc*Gp;   Lpi = Gpi*Gp;
T   = feedback(L,1);  S = 1/(1+L);   Tpi = feedback(Lpi,1);
out('\nFINAL  Gc(s) = Kp/(1+s/p), Kp = %.1f V/rad, p = %.2f rad/s  (Kv = Kp*Kdc = %.4f 1/s = 1/tp)\n', P.Kp_eff, P.p_eff, P.Kp_eff*Kdc);
[Gm,Pm,Wg,Wc] = margin(L);
out('  crossover wc = %.3f rad/s, PM = %.2f deg, GM = %.2f (%.2f dB) at %.2f rad/s\n', Wc, Pm, Gm, 20*log10(Gm), Wg);
bw = bandwidth(T);  out('  closed-loop -3 dB bandwidth = %.3f rad/s\n', bw);
Ms = norm(S, inf); Mt = norm(T, inf);
out('  Ms = %.3f (%.2f dB), Mt = %.3f\n', Ms, 20*log10(Ms), Mt);
si = stepinfo(T);  out('  step: overshoot %.2f %%, rise %.3f s, settling(2%%) %.3f s\n', si.Overshoot, si.RiseTime, si.SettlingTime);
cl = pole(T); out('  closed-loop poles: %s\n', mat2str(cl, 4));
out('  |T(j wn_plant)| = %.3e (%.1f dB) -> resonance attenuated\n', abs(freqresp(T, wn)), 20*log10(abs(freqresp(T, wn))));
% noise-equivalent bandwidth and 20 Hz sampled-noise reduction
w = logspace(-3, 3, 4000)';  Tm = abs(squeeze(freqresp(T, w)));
Bn = trapz(w, Tm.^2);  nred = sqrt(Bn/(pi/P.dt_meas));
out('  noise-equivalent BW = %.3f rad/s -> white 20 Hz reference noise reduced to %.3f x (%.1f dB)\n', Bn, nred, 20*log10(nred));
out('  Kv (ramp) = %.4f 1/s -> steady ramp lag = %.4f s (target tp = %.2f s)\n', P.Kp_eff*Kdc, 1/(P.Kp_eff*Kdc), P.tp);

% group delay of T at low frequency vs tp
wl = logspace(-2, 1, 600)';  ph = unwrap(angle(squeeze(freqresp(T, wl))));
gd = -gradient(ph, wl);
out('  group delay of T: %.4f s at 0.01 rad/s, %.4f s at 0.5 rad/s, %.4f s at 1 rad/s\n', gd(1), interp1(wl, gd, 0.5), interp1(wl, gd, 1));
wcmp = [0.1 0.25 0.5 1.0];
for q = 1:numel(wcmp)
    Tj = squeeze(freqresp(T, wcmp(q)));  dly = exp(-1j*wcmp(q)*P.tp);
    out('  w=%.2f: |T - e^{-j w tp}| = %.4f   (|T| = %.4f, phase err = %.2f deg)\n', wcmp(q), abs(Tj - dly), abs(Tj), rad2deg(angle(Tj/dly)));
end

out('\nSUGGESTED PI  Gc = 1800(1+1/s): ');
[Gm2,Pm2,Wg2,Wc2] = margin(Lpi);
out('wc = %.3f rad/s, PM = %.2f deg, GM = %.2f dB at %.2f rad/s; type-2 -> ramp lag = 0 (mount LEADS the object by tp)\n', Wc2, Pm2, 20*log10(Gm2), Wg2);
sip = stepinfo(Tpi); out('  step: overshoot %.2f %%, settling %.3f s\n', sip.Overshoot, sip.SettlingTime);
gdpi = -gradient(unwrap(angle(squeeze(freqresp(Tpi, wl)))), wl);
out('  group delay of T_pi at 0.1 rad/s = %.4f s, at 0.5 rad/s = %.4f s\n', interp1(wl,gdpi,0.1), interp1(wl,gdpi,0.5));

%% 3. Bode of L with margins ----------------------------------------------------------
f1 = figure('Position',[100 100 560 460],'Visible','off');
w = logspace(-1.5, 2.5, 800);
[Mg, Ph] = bode(L, w);  Mg = squeeze(Mg);  Ph = squeeze(Ph);
[Mg2, Ph2] = bode(Lpi, w); Mg2 = squeeze(Mg2); Ph2 = squeeze(Ph2);
ax1 = subplot(2,1,1);
semilogx(w, 20*log10(Mg), 'b', 'LineWidth', 1.6); hold on;
semilogx(w, 20*log10(Mg2), 'Color', [0.85 0.33 0.1], 'LineStyle', '--', 'LineWidth', 1.2);
yline(0, 'k:'); xline(Wc, 'k:'); xline(Wg, 'k:');
plot(Wg, -20*log10(Gm), 'rv', 'MarkerFaceColor','r');
text(Wg*1.1, -20*log10(Gm)+7, sprintf('GM = %.1f dB @ %.0f rad/s', 20*log10(Gm), Wg), 'FontSize', 8);
text(0.04, 12, sprintf('\\omega_c = %.2f rad/s', Wc), 'FontSize', 8);
ylabel('|L(j\omega)| [dB]'); grid on; legend('L = G_c G_p (final lag-matched)', 'PI (K_p=1800, T_i=1) [rejected]', 'Location','southwest','FontSize',7);
ylim([-80 50]); title('Open-loop Bode diagram with stability margins');
ax2 = subplot(2,1,2);
semilogx(w, Ph, 'b', 'LineWidth', 1.6); hold on;
semilogx(w, Ph2, 'Color', [0.85 0.33 0.1], 'LineStyle', '--', 'LineWidth', 1.2);
yline(-180, 'k:'); xline(Wc, 'k:'); xline(Wg, 'k:');
plot(Wc, Pm-180, 'rv', 'MarkerFaceColor','r');
text(Wc*1.2, Pm-180-22, sprintf('PM = %.1f^\\circ', Pm), 'FontSize', 8);
ylabel('\angle L [deg]'); xlabel('\omega [rad/s]'); grid on; ylim([-300 -60]);
linkaxes([ax1 ax2], 'x'); xlim([0.03 300]);
exportgraphics(f1, fullfile(figdir,'fig_bode.png'), 'Resolution', 300); exportgraphics(f1, fullfile(figdir,'fig_bode.pdf'));

%% 4. Closed-loop step, lag-matching, sensitivity -------------------------------------------
f2 = figure('Position',[100 100 900 300],'Visible','off');
subplot(1,3,1);
tt = linspace(0, 4, 800);
plot(tt, step(T, tt), 'b', 'LineWidth', 1.6); hold on; plot(tt, step(Tpi, tt), '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
yline(1, 'k:'); grid on; xlabel('t [s]'); ylabel('\theta/\theta_{ref}');
title('Unit step (linear model)'); legend('final','PI','Location','southeast','FontSize',7);
subplot(1,3,2);
semilogx(wl, gd, 'b', 'LineWidth', 1.6); hold on; semilogx(wl, gdpi, '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
yline(P.tp, 'k:'); grid on; xlabel('\omega [rad/s]'); ylabel('group delay of T [s]'); ylim([-0.1 0.8]);
title(sprintf('Loop lag vs t_p = %.1f s', P.tp));
subplot(1,3,3);
[Ms_, ~] = deal(0); %#ok<ASGLU>
w = logspace(-2, 2.5, 600); Sm_ = squeeze(abs(freqresp(S, w))); Tm_ = squeeze(abs(freqresp(T, w)));
semilogx(w, 20*log10(Tm_), 'b', 'LineWidth', 1.6); hold on; semilogx(w, 20*log10(Sm_), 'g', 'LineWidth', 1.2);
xline(wn, 'k:'); grid on; ylim([-80 10]); xlabel('\omega [rad/s]'); ylabel('[dB]'); legend('|T|','|S|','Location','southwest','FontSize',7);
title('Sensitivity and noise roll-off');
exportgraphics(f2, fullfile(figdir,'fig_closedloop.png'), 'Resolution', 300); exportgraphics(f2, fullfile(figdir,'fig_closedloop.pdf'));

%% 5. Tracking of a manoeuvring angle: lag-matched loop vs delayed reference ------------------
% reference: theta(t) with angular acceleration alpha; loop input = theta(t + tp); compare output with theta(t)
alpha = 0.025;                                    % rad/s^2  (= 20 m/s^2 at 800 m)
tt = (0:0.001:8)'; th0 = 0.3; om0 = 0.1;
thf = @(t) th0 + om0*t + 0.5*alpha*t.^2;
u  = thf(tt + P.tp);
y  = lsim(T, u, tt);
ey = y - thf(tt);
out('\nTracking a parabolic angle (alpha = %.3f rad/s^2 = 20 m/s^2 at 800 m): steady error y(t)-theta(t) = %.3e rad (%.3f m at 800 m, %.3f m at 2000 m)\n', ...
    alpha, ey(end), abs(ey(end))*800, abs(ey(end))*2000);
yp = lsim(Tpi, u, tt); eyp = yp - thf(tt);
out('   same input with PI: error = %.3e rad (%.1f m at 800 m)  [mount leads]\n', eyp(end), abs(eyp(end))*800);
out('   first-order-lag-only (no p): residual alpha*tp^2/2 = %.3e rad\n', alpha*P.tp^2/2);

%% 6. Saturation + anti-windup demonstration (nonlinear plant, 0.6 rad reference step) --------------------
dt = 1e-3; Tn = 8; nn = round(Tn/dt);
Xs = zeros(3,1); xf_aw = 0; xf_no = 0; Xn = zeros(3,1);
rec = zeros(nn, 6);
thref = 0.6;                                      % 34.4 deg step
for k = 1:nn
    t = (k-1)*dt;
    a = 1 - exp(-P.p_eff*dt);
    % with anti-windup (clamped state == output)
    e1 = thref - Xs(1);  xf_aw = xf_aw + a*(P.Kp_eff*e1 - xf_aw);  V1 = min(max(xf_aw,-P.Vmax),P.Vmax); xf_aw = V1;
    % without anti-windup (state evolves unsaturated, only the output is clipped)
    e2 = thref - Xn(1);  xf_no = xf_no + a*(P.Kp_eff*e2 - xf_no);  V2 = min(max(xf_no,-P.Vmax),P.Vmax);
    rec(k,:) = [t, Xs(1), V1, Xn(1), V2, xf_no];
    Xs = rk4(Xs, V1, dt); Xn = rk4(Xn, V2, dt);
end
f3 = figure('Position',[100 100 700 300],'Visible','off');
subplot(1,2,1); plot(rec(:,1), rad2deg(rec(:,2)), 'b', 'LineWidth', 1.5); hold on; plot(rec(:,1), rad2deg(rec(:,4)), '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
yline(rad2deg(thref), 'k:'); grid on; xlabel('t [s]'); ylabel('\theta_{turret} [deg]'); legend('anti-windup (state clamp)', 'no anti-windup', 'Location','southeast','FontSize',7);
title('34.4^\circ reference step, nonlinear plant');
subplot(1,2,2); plot(rec(:,1), rec(:,3), 'b', 'LineWidth', 1.5); hold on; plot(rec(:,1), rec(:,5), '--', 'Color', [0.85 0.33 0.1], 'LineWidth', 1.2);
yline(P.Vmax, 'k:'); yline(-P.Vmax, 'k:'); grid on; xlabel('t [s]'); ylabel('V_m [V]'); title('Motor voltage (\pm240 V saturation)');
exportgraphics(f3, fullfile(figdir,'fig_antiwindup.png'), 'Resolution', 300); exportgraphics(f3, fullfile(figdir,'fig_antiwindup.pdf'));
os_aw = (max(rec(:,2))/thref - 1)*100;  os_no = (max(rec(:,4))/thref - 1)*100;
t98 = @(c) rec(find(abs(rec(:,c)-thref) > 0.02*thref, 1, 'last')+1, 1);
out('\nSaturated 0.6 rad step: overshoot %.2f %% with anti-windup (2%% settle %.2f s) vs %.2f %% without (2%% settle %.2f s)\n', os_aw, t98(2), os_no, t98(4));
out('   slew-limited rate at saturation = %.4f rad/s (theory %.4f)\n', max(rec(1000:2500,2)-rec(999:2499,2))/dt, P.Vmax*Kdc);
fclose(fid);
save(fullfile(resdir,'design_numbers.mat'),'Wc','Pm','Gm','Wg','bw','Ms','Mt','si','Bn','nred','wn','zeta');
% Wc Pm Gm Wg bw Ms Mt si Bn nred wn zeta
disp('design_analysis done.');

function X = rk4(X, V, dt)
f = @(X) local_f(V, X);
k1 = f(X); k2 = f(X+dt/2*k1); k3 = f(X+dt/2*k2); k4 = f(X+dt*k3);
X = X + dt/6*(k1+2*k2+2*k3+k4);
end
function d = local_f(V, X)
[a,b,c] = turret_plant(0, V, X(1), X(2), X(3)); d = [a;b;c];
end

function fprintf_both(fid, fmt, varargin)
fprintf(fid, fmt, varargin{:}); fprintf(fmt, varargin{:});
end
