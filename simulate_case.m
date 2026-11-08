function L = simulate_case(scenario, ovr, ctrl_dt)
%SIMULATE_CASE  Closed-loop simulation of one case over the FULL 0 <= t <= 10 s.
%   L = simulate_case(scenario)        scenario: 'A','B','C' or struct('seed',n[,'loss',..])
%   L = simulate_case(scenario, ovr)   ovr: struct of turret_params overrides (ablations only)
%   L = simulate_case(scenario, ovr, ctrl_dt)  controller call period (default 1 ms;
%       must be a multiple of 1 ms; a measurement that arrives between two calls is
%       latched and delivered at the next call)
%
%   Fixed-step RK4 (1 ms) of the plant [theta, theta_dot, i]; the controller is
%   called every step (zero-order hold on V_m).  Pulse logic (driver side):
%   a pulse is issued on a 0->1 edge of fire_cmd if >= 1.0 s since the last
%   issued pulse; premature edges are ignored.  A pulse succeeds when the true
%   perpendicular distance object-to-beam is <= R_eng = 2 m with the object on
%   the forward side.
clear turret_controller target_motion_generator turret_params turret_plant
if nargin > 1 && ~isempty(ovr), turret_params(ovr); end
P = turret_params();
dt = 1e-3;  Tend = 10;  n = round(Tend/dt);
if nargin < 3 || isempty(ctrl_dt), ctrl_dt = dt; end
mstep = round(ctrl_dt/dt);  rm_l = NaN; thm_l = NaN; pend = false;
X = [0; 0; 0];                                    % theta, theta_dot, i  (all zero at t = 0)

N = n + 1;
f = {'t','x','y','vx','vy','ax','ay','th','thd','cur','V','Vun','xh','yh','vxh','vyh','axh','ayh', ...
     'sx','sy','svx','svy','thref','dperp','dperp_est','xpx','xpy','loss','sigperp'};
for k = 1:numel(f), L.(f{k}) = nan(N,1); end
pulses = zeros(0,4);                              % [t, d_perp_true, success, range]
t_last = -inf; fire_prev = 0;
RK = @(t,V,X) local_deriv(t,V,X);

for k = 0:n
    t = k*dt;  j = k+1;
    [rm, thm, valid, tr] = target_motion_generator(t, scenario);
    if valid, rm_l = rm; thm_l = thm; pend = true; end          % latch new measurement
    if mod(k, mstep) == 0                                        % controller sample instant
        [V, fire, dbg] = turret_controller(t, rm_l, thm_l, pend, X(1), X(2));
        pend = false;
    else
        fire = 0;                                                % fire_cmd is a one-sample signal
    end
    % ---- driver-side pulse handling
    if fire && ~fire_prev && (t - t_last) >= P.t_refr - 1e-9
        c = cos(X(1)); s = sin(X(1));
        dtrue = abs(tr.x*s - tr.y*c);  fw = tr.x*c + tr.y*s;
        ok = (dtrue <= P.R_eng) && (fw > 0);
        pulses(end+1,:) = [t, dtrue, ok, hypot(tr.x,tr.y)]; %#ok<AGROW>
        t_last = t;
    end
    fire_prev = fire;
    % ---- log
    L.t(j)=t; L.x(j)=tr.x; L.y(j)=tr.y; L.vx(j)=tr.vx; L.vy(j)=tr.vy; L.ax(j)=tr.ax; L.ay(j)=tr.ay;
    L.th(j)=X(1); L.thd(j)=X(2); L.cur(j)=X(3); L.V(j)=V; L.Vun(j)=dbg.V_unsat;
    L.xh(j)=dbg.x(1); L.yh(j)=dbg.x(2); L.vxh(j)=dbg.x(3); L.vyh(j)=dbg.x(4); L.axh(j)=dbg.x(5); L.ayh(j)=dbg.x(6);
    L.sx(j)=dbg.sigma(1); L.sy(j)=dbg.sigma(2); L.svx(j)=dbg.sigma(3); L.svy(j)=dbg.sigma(4);
    L.thref(j)=dbg.th_ref; L.dperp_est(j)=dbg.dperp; L.xpx(j)=dbg.xp(1); L.xpy(j)=dbg.xp(2);
    L.sigperp(j)=dbg.sig_perp; L.loss(j)=tr.in_loss;
    c = cos(X(1)); s = sin(X(1));
    L.dperp(j) = abs(tr.x*s - tr.y*c);
    if (tr.x*c + tr.y*s) <= 0, L.dperp(j) = hypot(tr.x,tr.y); end   % behind the mount: no valid beam
    if k == n, break; end
    % ---- RK4 plant step with V held
    k1 = RK(t, V, X);  k2 = RK(t+dt/2, V, X+dt/2*k1);
    k3 = RK(t+dt/2, V, X+dt/2*k2);  k4 = RK(t+dt, V, X+dt*k3);
    X = X + dt/6*(k1 + 2*k2 + 2*k3 + k4);
    if X(1) > P.th_max, X(1) = P.th_max; X(2) = min(X(2),0); end
    if X(1) < P.th_min, X(1) = P.th_min; X(2) = max(X(2),0); end
end
L.pulses = pulses;  L.scenario = scenario;  L.P = P;
L.theta_los = atan2(L.y, L.x);  L.r = hypot(L.x, L.y);
L.metrics = local_metrics(L);
end

function dX = local_deriv(t, V, X)
[a, b, c] = turret_plant(t, V, X(1), X(2), X(3));
dX = [a; b; c];
end

function M = local_metrics(L)
P = L.P;  dt = L.t(2) - L.t(1);  tp = P.tp;
M.rms_pos  = sqrt(mean((L.xh-L.x).^2 + (L.yh-L.y).^2));                       % m, full 0..10 s
M.rms_vel  = sqrt(mean((L.vxh-L.vx).^2 + (L.vyh-L.vy).^2));                   % m/s
% prediction error  |x_pred(t) - x_true(t+tp)|  for t <= 10 - tp
np = round(tp/dt);  m = 1:(numel(L.t)-np);
pred_err = sqrt((L.xpx(m)-L.x(m+np)).^2 + (L.xpy(m)-L.y(m+np)).^2);
M.rms_pred = sqrt(mean(pred_err.^2));
% prediction error restricted to the coast-stressed window: during a loss and the
% 0.5 s after it (the filter is extrapolating open-loop along a possibly-curved path)
coast = false(size(L.t));
li = find(diff([0; L.loss>0]) ~= 0);                 % loss edges
lw = L.loss > 0;
for q = 1:numel(L.t)
    if lw(q), coast(q) = true; end
end
coast_ext = coast;
for q = 2:numel(L.t)                                  % extend 0.5 s past each loss
    if coast(q-1) && ~coast(q)
        coast_ext(q:min(q+round(0.5/dt),numel(L.t))) = true;
    end
end
cm = coast_ext(m);
if any(cm), M.rms_pred_loss = sqrt(mean(pred_err(cm).^2)); else, M.rms_pred_loss = NaN; end
% position estimation error during loss windows only
if any(lw), M.rms_pos_loss = sqrt(mean((L.xh(lw)-L.x(lw)).^2 + (L.yh(lw)-L.y(lw)).^2));
else, M.rms_pos_loss = NaN; end
% true LOS angular rate (from truth) -- used for the high-rate pointing metric
r2 = L.x.^2 + L.y.^2;  thd_los = (L.x.*L.vy - L.y.*L.vx)./max(r2,1e-9);
% elevation errors
e_ref = L.thref - L.th;                         % loop tracking error
e_los = L.theta_los - L.th;                     % pointing error w.r.t. the object's actual LOS
M.rms_eref_full = sqrt(mean(e_ref.^2));         % rad (includes initial slew)
M.rms_elos_full = sqrt(mean(e_los.^2));
% acquisition time: last time the true miss distance was > 10 m
idx = find(L.dperp > 10, 1, 'last');
if isempty(idx), tacq = 0; elseif idx == numel(L.t), tacq = NaN; else, tacq = L.t(idx+1); end
M.t_acq = tacq;
if isfinite(tacq)
    q = L.t >= tacq;
    M.rms_eref_acq = sqrt(mean(e_ref(q).^2));
    M.rms_elos_acq = sqrt(mean(e_los(q).^2));
    M.max_elos_acq = max(abs(e_los(q)));                        % peak pointing error post-acq (rad)
    % pointing error during HIGH angular-rate windows (where lag-match breaks down)
    hr = q & (abs(thd_los) > 0.10);                            % |LOS rate| > 0.10 rad/s
    if any(hr), M.rms_elos_hirate = sqrt(mean(e_los(hr).^2)); else, M.rms_elos_hirate = NaN; end
    M.max_miss_acq = max(L.dperp(q));               % NOTE: by construction this sits near the 10 m
                                                      % crossing threshold itself (it is the boundary
                                                      % point of the "acquired" definition) -- use
                                                      % max_miss_acq_buf / rms_miss_acq for the actual
                                                      % steady-state accuracy reported in the paper.
    M.rms_miss_acq = sqrt(mean(L.dperp(q).^2));
    qb = L.t >= (tacq + 0.75);                       % buffered window: excludes the crossing transient
    if any(qb)
        M.max_miss_acq_buf = max(L.dperp(qb));
        M.rms_miss_acq_buf = sqrt(mean(L.dperp(qb).^2));
    else
        M.max_miss_acq_buf = NaN; M.rms_miss_acq_buf = NaN;
    end
else
    M.rms_eref_acq = NaN; M.rms_elos_acq = NaN; M.max_miss_acq = NaN; M.rms_miss_acq = NaN;
    M.max_miss_acq_buf = NaN; M.rms_miss_acq_buf = NaN;
    M.max_elos_acq = NaN; M.rms_elos_hirate = NaN;
end
M.n_pulses = size(L.pulses,1);
M.n_success = sum(L.pulses(:,3) == 1);
M.success = M.n_success > 0;
if M.success, M.t_first_hit = L.pulses(find(L.pulses(:,3)==1,1),1); else, M.t_first_hit = NaN; end
if M.n_pulses > 0, M.max_miss_pulse = max(L.pulses(:,2)); M.first_pulse_miss = L.pulses(1,2);
else, M.max_miss_pulse = NaN; M.first_pulse_miss = NaN; end
M.Vsat_frac = mean(abs(L.Vun) > P.Vmax);        % fraction of time the command is saturated
M.max_theta = max(L.th); M.min_theta = min(L.th);
end
