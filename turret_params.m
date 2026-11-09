function p = turret_params(ovr)
%TURRET_PARAMS  Single source of truth for every constant of the design.
%   p = turret_params()   returns the persistent parameter struct.
%   turret_params(ovr)    overrides fields of the struct (used ONLY by the
%                         ablation studies in run_headtohead.m).
%   clear turret_params   restores the defaults.
%
%   The default values are THE finalized design.  They are never changed between
%   scenarios or evaluation cases: every alternative reported in the paper is a
%   separate, frozen configuration reached through an override.
%
%   Notation follows the brief (Assignment 2, ME319 Autumn 2026, Sec. 2):
%     turret elevation theta_t, motor voltage V_m, motor current i,
%     range r, elevation angle theta, measurement-loss interval t_bo <= 0.50 s
%     starting no earlier than t_boo >= 1.00 s, engagement radius R_eng = 2 m,
%     recharge time t_recharge = 1.0 s, prediction horizon Tp = 0.50 s.
%
%   Units: SI (m, s, rad, V, A, N m).
persistent P
if isempty(P), P = local_defaults(); end
if nargin > 0 && isstruct(ovr)
    f = fieldnames(ovr);
    for k = 1:numel(f), P.(f{k}) = ovr.(f{k}); end
    P = local_derive(P);
end
p = P;
end

% =========================================================================
function P = local_defaults()
% ---- elevation-drive plant (given in the brief, Sec. 2) -------------------
P.J_turret = 1200;      % kg m^2   turret + cradle inertia about the trunnion
P.J_m      = 0.20;      % kg m^2   motor rotor inertia
P.N        = 100;       % -        gear ratio (motor speed / turret speed)
P.b        = 1500;      % N m s/rad turret-side viscous rotational damping
P.Kt       = 10;        % N m/A    motor torque constant
P.Ke       = 10;        % V s/rad  back-emf constant
P.Rm       = 1.0;       % ohm      armature resistance
P.Lm       = 0.10;      % H        armature inductance
P.th_min   = deg2rad(-10);   % rad  mechanical elevation limit
P.th_max   = deg2rad(85);    % rad  mechanical elevation limit
P.Vmax     = 240;       % V        actuator (voltage) saturation, assumed
% ---- radar / operating envelope (given in the brief, Secs. 2-3) -----------
P.dt_meas  = 0.05;      % s       20 Hz measurement interval
P.sig_r    = 1.0;       % m       range noise 1-sigma
P.sig_th   = 0.02*pi/180;   % rad elevation-angle noise 1-sigma (= 0.02 deg)
P.r_min = 800; P.r_max = 2000;   % m  operating range
P.v_max = 150; P.a_max = 20;     % m/s, m/s^2  target speed / acceleration limits
P.t_bo_max  = 0.50;     % s       maximum measurement-loss duration
P.t_boo_min = 1.00;     % s       earliest measurement-loss onset
P.T_end     = 10;       % s       standardized engagement duration
P.dt_sim    = 1e-3;     % s       RK4 integration step
% ---- engagement (laser) logic --------------------------------------------
P.R_eng    = 2.0;       % m       engagement radius of a pulse
P.t_recharge = 1.0;     % s       minimum interval between engagement pulses
P.fire_thr = 1.5;       % m       fire when the estimated perpendicular miss <= this
P.fire_sig = 1.0;       % m       ... and the 1-sigma cross-beam uncertainty <= this
% ---- target-state estimator / predictor ----------------------------------
P.q_jerk   = 200;       % m^2/s^5 jerk spectral density of the CA white-jerk model
P.sig_v0   = 100;       % m/s     initial velocity 1-sigma (deliberately large)
P.sig_a0   = 15;        % m/s^2   initial acceleration 1-sigma (deliberately large)
P.tp       = 0.50;      % s       prediction horizon Tp (FIXED by the brief)
P.pred_model = 'CA';    % 'CA' (default) | 'CV' (ablation only)
% ---- estimator selector (default = the finalized single-model CA filter) ---
P.est      = 'CA';      % 'CA' (final) | 'adaptQ' (NIS-adaptive Q) | 'IMM' (CV/CA/CT)
%   adaptive-Q: declare a manoeuvre when the normalized innovation squared (NIS,
%   chi-square, 2 dof) exceeds aq_gate and transiently inflate the jerk PSD.
P.aq_gate  = 9.21;      % chi2inv(0.99,2): NIS gate for "manoeuvre declared"
P.aq_boost = 50;        % jerk-PSD multiplier applied on detection
P.aq_relax = 0.6;       % geometric decay of the multiplier toward 1 per measurement
%   IMM: three models in a common [x y vx vy ax ay] frame.
P.imm_q_cv    = 3;      % m^2/s^3   CV (white-acceleration) spectral density
P.imm_q_ca    = 200;    % m^2/s^5   CA jerk PSD (matches the single-model baseline)
P.imm_q_ct    = 150;    % m^2/s^5   CT base jerk PSD
P.imm_q_omega = 6.0;    % (rad/s^2)^2 extra process variance on the turn-rate state
P.imm_TPM     = [0.95 0.025 0.025;    % CV
                 0.025 0.95 0.025;    % CA
                 0.025 0.025 0.95];   % CT   (row = from, col = to)
% ---- feedback compensator -------------------------------------------------
%   FINALIZED controller: 2-DOF LTI = lag-matched feedback + model-inverse
%   feedforward.  The feedback is the conventional LTI compensator required by
%   the brief; the feedforward is an LTI, analytically DERIVED (not tuned) path
%   that holds the closed-loop group delay at Tp under angular acceleration.
%   Ablations: 'lag' (feedback only), 'pi' (the brief's suggested PI).
P.ctrl      = 'lagff';  % 'lagff' = lag feedback + 2-DOF feedforward (FINAL) | 'lag' | 'pi'
P.tp_design = 0.50;     % s       group delay the loop is DESIGNED to have (= Tp)
P.Kp_pi = 1800; P.Ti_pi = 1.0;   % suggested PI starting values (ablation only)
P.Kp_lag = []; P.p_lag = [];     % derived below unless explicitly overridden
P.tp_ff  = 0.0;         % s       aim horizon used with the feedforward: the FF path
                        %         supplies the Tp lead itself, so the feedback sees the
                        %         current estimate; 'lag' alone aims at Tp
% ---- acquisition front-end ------------------------------------------------
%   FINALIZED: bounded dual-mode acquisition (rate-limited slew, then a bumpless
%   hand-off to the LTI tracker).  The engagement-phase controller is the LTI
%   2-DOF design above.  Ablation: 'linear' (LTI tracker throughout).
P.acq      = 'dualmode';% 'dualmode' (FINAL) | 'linear' (ablation)
P.acq_band = deg2rad(2.5);   % rad  slew -> fine tracking hand-off band.
                             %   Chosen from the geometry of the saturated-rate unwind:
                             %   the mount leaves saturation at ~0.24 rad/s, and the
                             %   compensator's 1/p = 0.25 s relaxation covers ~3.4 deg,
                             %   so a band below ~1.5 deg re-enters saturation and the
                             %   capture oscillates; 2.5 deg keeps the hand-off clean
                             %   while leaving the slew phase as short as possible.
P.acq_Kp   = 20000;     % V/rad   slew-phase gain: always saturates -> rate-limited slew
P = local_derive(P);
end

% =========================================================================
function P = local_derive(P)
P.J_eq = P.J_turret + P.N^2*P.J_m;                              % 3200 kg m^2
% Low-frequency velocity gain of the plant, d(theta)/dt = Kdc*V   [rad/s per V].
% Obtained from Gp(s) with s -> 0 after cancelling the integrator:
P.Kdc  = P.N*P.Kt / (P.N^2*P.Kt*P.Ke + P.Rm*P.b);               % 9.985e-4
% Lag-matched feedback design (paper Sec. IV-A):  Kv = Kp*Kdc = 1/Tp,  p = 2/Tp.
P.Kp_eff = P.Kp_lag; P.p_eff = P.p_lag;
if isempty(P.Kp_eff), P.Kp_eff = 1/(P.tp_design*P.Kdc); end     % 2003 V/rad
if isempty(P.p_eff),  P.p_eff  = 2/P.tp_design;           end    % 4 rad/s
% 2-DOF feedforward gains (model inverse, DERIVED not tuned):
%   V_ff = thetadot_ref/Kdc  drives the turret at the reference's own angular
%   rate (back-emf term); the acceleration term supplies the armature-resistance
%   voltage needed to source the current for J_eq*thetaddot of torque.
P.kff_v = 1/P.Kdc;                       % V per rad/s   (rate feedforward)
P.kff_a = P.Rm*P.J_eq/(P.N*P.Kt);        % V per rad/s^2 (acceleration feedforward)
end
