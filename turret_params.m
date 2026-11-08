function p = turret_params(ovr)
%TURRET_PARAMS  Single source of truth for every design constant.
%   p = turret_params()      returns the (persistent) parameter struct.
%   turret_params(ovr)       overrides fields of the struct (used ONLY by the
%                            ablation studies in run_final_evaluation.m).
%   clear turret_params      restores the defaults.
%
%   The default values below are THE finalized design.  They are not changed
%   between scenarios or evaluation cases.
persistent P
if isempty(P), P = local_defaults(); end
if nargin > 0 && isstruct(ovr)
    f = fieldnames(ovr);
    for k = 1:numel(f), P.(f{k}) = ovr.(f{k}); end
    P = local_derive(P);
end
p = P;
end

% -------------------------------------------------------------------------
function P = local_defaults()
% ---- plant (Sec. 3 of the assignment) -----------------------------------
P.J_turret = 1200;      % kg m^2   mount + cradle inertia
P.J_m      = 0.20;      % kg m^2   rotor inertia
P.N        = 100;       % -        gear ratio (motor speed / mount speed)
P.b        = 1500;      % N m s/rad mount-side viscous damping
P.Kt       = 10;        % N m/A
P.Ke       = 10;        % V s/rad
P.Rm       = 1.0;       % ohm
P.Lm       = 0.10;      % H
P.th_min   = deg2rad(-10);
P.th_max   = deg2rad(85);
P.Vmax     = 240;       % V   actuator saturation
% ---- sensor / envelope ----------------------------------------------------
P.dt_meas  = 0.05;                 % s   (20 Hz)
P.sig_r    = 1.0;                  % m
P.sig_th   = 0.02*pi/180;          % rad (3.4907e-4)
P.r_min = 800; P.r_max = 2000; P.v_max = 150; P.a_max = 20;
% ---- pulse logic -----------------------------------------------------------
P.R_eng    = 2.0;       % m     success radius
P.t_refr   = 1.0;       % s     refractory period
P.fire_thr = 1.5;       % m     fire when estimated miss <= this (margin 0.5 m)
P.fire_sig = 1.0;       % m     ... and 1-sigma cross-beam position uncertainty <= this
% ---- estimator / predictor -------------------------------------------------
P.q_jerk   = 200;       % m^2/s^5  jerk spectral density (CA white-jerk model)
P.sig_v0   = 100;       % m/s   initial velocity 1-sigma (large)
P.sig_a0   = 15;        % m/s^2 initial acceleration 1-sigma (large)
P.tp       = 0.5;       % s     prediction horizon (FIXED by the assignment)
P.pred_model = 'CA';    % 'CA' (default) or 'CV' (ablation only)
% ---- compensator -------------------------------------------------------------
P.ctrl      = 'lag';    % 'lag' = Kp/(1+s/p)  (final);  'pi' = suggested PI (ablation)
P.tp_design = 0.5;      % s     lag the loop is DESIGNED to have (= tp)
P.Kp_pi = 1800; P.Ti_pi = 1.0;     % suggested PI starting values (ablation only)
P.Kp_lag = []; P.p_lag = [];       % derived below unless overridden
P = local_derive(P);
end

function P = local_derive(P)
P.J_eq = P.J_turret + P.N^2*P.J_m;                              % 3200 kg m^2
% low-frequency velocity gain of the plant:  dtheta/dt = Kdc * V   (rad/s per V)
P.Kdc  = P.N*P.Kt / (P.N^2*P.Kt*P.Ke + P.Rm*P.b);               % 9.985e-4
% Lag-matched design (see report Sec. IV):  Kv = Kp*Kdc = 1/tp,  p = 2/tp
P.Kp_eff = P.Kp_lag; P.p_eff = P.p_lag;
if isempty(P.Kp_eff), P.Kp_eff = 1/(P.tp_design*P.Kdc); end    % 2003 V/rad
if isempty(P.p_eff),  P.p_eff  = 2/P.tp_design;          end    % 4 rad/s
end
