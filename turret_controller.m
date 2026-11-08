function [V_m, fire_cmd, dbg] = turret_controller(t, r_meas, theta_meas, meas_valid, theta_turret, theta_dot_turret) %#ok<INUSD>
%TURRET_CONTROLLER  Predictive tracking controller for the elevation mount.
%
%   Data flow (Sec. 7):  polar measurement -> target-state estimator (Cartesian)
%     -> tp = 0.5 s prediction -> theta_ref = atan2(y_pred, x_pred) (clamped)
%     -> LTI compensator (+/- 2-DOF feedforward, saturation, anti-windup)
%     -> V_m, plus single-sample fire_cmd logic with a 1 s refractory period.
%
%   The estimator, compensator and acquisition mode are selected by turret_params:
%     P.est  = 'CA' (baseline) | 'adaptQ' (NIS-adaptive jerk PSD) | 'IMM' (CV/CA/CT)
%     P.ctrl = 'lag' (baseline) | 'pi' (ablation) | 'lagff' (lag + 2-DOF feedforward)
%     P.acq  = 'linear' (baseline) | 'dualmode' (near-time-optimal slew + capture)
%   With the default parameters the behaviour is IDENTICAL to the baseline design.
%
%   Persistent state (cleared with `clear turret_controller`; also re-initialised
%   automatically if time ever decreases): estimator state/covariance (and, for
%   IMM, the sub-filter bank + model probabilities), compensator state, last
%   pulse time, previous fire_cmd.
%
%   The controller may be called at ANY rate (all discretisations use the actual
%   elapsed time t - t_prev).  theta_dot_turret is part of the required interface
%   but is not used by the final (single-loop) design.
persistent S
P = turret_params();
if isempty(S) || t < S.t_prev - 1e-9
    S = local_new_state(t, P);
end
dt = t - S.t_prev;

%% ---------------------------------------------------------------- estimator
if ~S.init
    if meas_valid && isfinite(r_meas) && isfinite(theta_meas)
        S = local_est_init(S, r_meas, theta_meas, P);
        S.init = true;
    end
else
    S = local_est_step(S, dt, meas_valid, r_meas, theta_meas, P);
end
S.t_prev = t;

%% --------------------------------------------- prediction -> elevation reference
use_ff = strcmp(P.ctrl, 'lagff');
if S.init
    xp = local_est_predict(S, P.tp, P);           % ALWAYS the tp=0.5 s prediction (for metric/pulse-free)
    h_aim = P.tp; if use_ff, h_aim = P.tp_ff; end  % aim horizon (FF supplies the lead otherwise from tp)
    if abs(h_aim - P.tp) < 1e-12, xa = xp; else, xa = local_est_predict(S, h_aim, P); end
    [th_ref, thref_dot, thref_ddot] = local_refkin(xa);
    th_ref_u = th_ref;
    th_ref = min(max(th_ref, P.th_min), P.th_max);        % mount limits
    if th_ref ~= th_ref_u, thref_dot = 0; thref_ddot = 0; end   % clamped -> no FF lead
else
    xp = [NaN;NaN;NaN;NaN;NaN;NaN]; th_ref = theta_turret; thref_dot = 0; thref_ddot = 0;
end

%% --------------------------------- 2-DOF feedforward voltage (0 unless ctrl='lagff')
V_ff = 0;
if use_ff && S.init
    V_ff = P.kff_v*thref_dot + P.kff_a*thref_ddot;
end

%% ----------------------------------------------- LTI compensator + saturation
e = th_ref - theta_turret;
slew = false;
if S.init && strcmp(P.acq,'dualmode') && abs(e) > P.acq_band
    slew = true;                                   % --- near-time-optimal max-rate slew phase
    V_unsat = P.acq_Kp*e + V_ff;                   % high gain -> saturates; FF keeps the handoff bumpless
    V_m = min(max(V_unsat, -P.Vmax), P.Vmax);
    S.xf = V_m - V_ff;                             % keep lag state consistent for a bumpless handoff
    S.xi = 0;
elseif ~S.init
    V_unsat = 0; V_m = 0;
elseif strcmp(P.ctrl, 'pi')                       % ablation: suggested PI
    xi_new = S.xi + e*dt;
    V_unsat = P.Kp_pi*(e + xi_new/P.Ti_pi);
    V_m = min(max(V_unsat, -P.Vmax), P.Vmax);
    if V_m ~= V_unsat                              % back-calculation anti-windup
        xi_new = P.Ti_pi*(V_m/P.Kp_pi - e);
    end
    S.xi = xi_new;
else                                              % 'lag' (baseline) or 'lagff' (lag + feedforward)
    a = 1 - exp(-P.p_eff*dt);                      % exact first-order discretisation
    xf = S.xf + a*(P.Kp_eff*e - S.xf);
    V_unsat = xf + V_ff;
    V_m = min(max(V_unsat, -P.Vmax), P.Vmax);
    if use_ff
        S.xf = V_m - V_ff;                         % back-calculation anti-windup around the FF term
    else
        S.xf = V_m;                                % baseline: state == output -> exact anti-windup
    end
end
% never drive into a mechanical stop
if (theta_turret >= P.th_max && V_m > 0) || (theta_turret <= P.th_min && V_m < 0)
    V_m = 0; S.xf = 0;
end

%% -------------------------------------------------------------- pulse logic
fire_cmd = 0; dperp = NaN; sig_perp = NaN; fwd = NaN;
if S.init
    c = cos(theta_turret); s = sin(theta_turret);
    px = S.x(1); py = S.x(2);
    fwd   = px*c + py*s;                           % forward projection
    dperp = abs(px*s - py*c);                      % perpendicular miss (estimated)
    u = [-s; c];  sig_perp = sqrt(u' * S.P(1:2,1:2) * u);
end
if S.fire_prev
    fire_cmd = 0; S.fire_prev = false;             % guarantee a clean rising edge
elseif S.init && fwd > 0 && dperp <= P.fire_thr && sig_perp <= P.fire_sig ...
        && (t - S.t_fire) >= P.t_refr - 1e-9
    fire_cmd = 1; S.fire_prev = true; S.t_fire = t;
end

if nargout > 2
    dbg.x = S.x; dbg.sigma = sqrt(max(diag(S.P), 0)); dbg.xp = xp(1:2);
    dbg.th_ref = th_ref; dbg.V_unsat = V_unsat; dbg.dperp = dperp; dbg.sig_perp = sig_perp;
    dbg.init = S.init; dbg.thref_dot = thref_dot; dbg.V_ff = V_ff; dbg.slew = slew;
    dbg.nis = S.nis; dbg.mu = S.mu;
end
end

% =========================================================================
% State container
% =========================================================================
function S = local_new_state(t, P)
S.init = false; S.x = zeros(6,1); S.P = eye(6);
S.t_prev = t; S.xf = 0; S.xi = 0;
S.t_fire = -inf; S.fire_prev = false;
S.qmult = 1; S.nis = 0;                            % adaptive-Q bookkeeping
S.mu = [0 0 0];                                    % IMM model probabilities (reported)
S.imm = [];                                        % IMM sub-filter bank
end

% =========================================================================
% Unified estimator interface:  init / step / predict
% =========================================================================
function S = local_est_init(S, r, th, P)
switch P.est
    case 'IMM'
        [x0, P0] = local_kf_init(r, th, P);
        mdl = {'cv','ca','ct'};
        for j = 1:3
            [xn, Pn] = local_from_common(mdl{j}, x0, P0);
            S.imm(j).x = xn; S.imm(j).P = Pn; S.imm(j).type = mdl{j};
        end
        S.mu = [1 1 1]/3;
        S.x = x0; S.P = P0;
    otherwise                                      % 'CA' (baseline) or 'adaptQ'
        [S.x, S.P] = local_kf_init(r, th, P);
        S.qmult = 1;
end
end

function S = local_est_step(S, dt, meas_valid, r, th, P)
has_meas = meas_valid && isfinite(r) && isfinite(th);
switch P.est
    case 'IMM'
        S = local_imm_step(S, dt, has_meas, r, th, P);
    case 'adaptQ'
        if dt > 0
            [S.x, S.P] = local_kf_predict(S.x, S.P, dt, P, P.q_jerk*S.qmult);
        end
        if has_meas
            [S.x, S.P, nis] = local_kf_update(S.x, S.P, r, th, P);
            S.nis = nis;
            if nis > P.aq_gate, S.qmult = P.aq_boost;          % manoeuvre declared
            else, S.qmult = 1 + (S.qmult - 1)*P.aq_relax; end  % relax toward baseline
        end
    otherwise                                      % 'CA' baseline (exact legacy path)
        if dt > 0, [S.x, S.P] = local_kf_predict(S.x, S.P, dt, P, P.q_jerk); end
        if has_meas, [S.x, S.P] = local_kf_update(S.x, S.P, r, th, P); end
end
end

function xp = local_est_predict(S, h, P)
% Combined 6-state prediction at horizon h (model-appropriate).
switch P.est
    case 'IMM'
        xp = zeros(6,1);
        for j = 1:3
            ty = S.imm(j).type;
            xn = local_model_fn(ty, S.imm(j).x, h);          % native prediction over h
            xc = local_to_common(ty, xn);                    % -> common 6-state
            xp = xp + S.mu(j) * xc;
        end
    otherwise
        if strcmp(P.pred_model,'CV')
            xp = S.x;  xp(1:2) = S.x(1:2) + S.x(3:4)*h;  xp(5:6) = 0;
        else
            xp = local_F(h) * S.x;                  % CA prediction
        end
end
end

% =========================================================================
% IMM (CV + CA in the common [x y vx vy ax ay] frame; CT as a native 5-state
% [x y vx vy omega] EKF with turn rate as a FILTERED state, mixed via common
% frame maps).  Reverts to the CA single model on non-turning motion and
% follows the arc on sustained coordinated turns.
% =========================================================================
function S = local_imm_step(S, dt, has_meas, r, th, P)
M = 3; mu = S.mu(:); TPM = P.imm_TPM;
% ---- 1. mixing probabilities
cbar = TPM' * mu;                                   % M x 1
mix = (TPM .* (mu * (1./max(cbar,1e-12))'));        % mix(i,j) = pi_ij mu_i / cbar_j
% ---- bring every model into the common 6-frame for mixing
XC = cell(M,1); PC = cell(M,1);
for i = 1:M, [XC{i}, PC{i}] = local_to_common(S.imm(i).type, S.imm(i).x, S.imm(i).P); end
% ---- 2. mixed initial conditions (common frame), then back to native
L = ones(M,1); nis_rep = 0;
for j = 1:M
    xj = zeros(6,1);
    for i = 1:M, xj = xj + mix(i,j)*XC{i}; end
    Pj = zeros(6);
    for i = 1:M, d = XC{i} - xj; Pj = Pj + mix(i,j)*(PC{i} + d*d'); end
    ty = S.imm(j).type;
    [xn, Pn] = local_from_common(ty, xj, Pj);       % native mixed init
    % ---- 3. model-matched predict + (optional) update, in NATIVE frame
    if dt > 0, [xn, Pn] = local_model_predict_n(ty, xn, Pn, dt, P); end
    if has_meas
        [xn, Pn, nisj, Lj] = local_update_n(ty, xn, Pn, r, th, P);
        L(j) = Lj;
        if strcmp(ty,'ca'), nis_rep = nisj; end     % report CA-model NIS for plotting
    end
    S.imm(j).x = xn; S.imm(j).P = Pn;
end
% ---- 4. model-probability update
if has_meas, mu = cbar .* L; else, mu = cbar; end   % coast: no likelihood information
mu = mu / max(sum(mu), 1e-300);
S.mu = mu(:)'; S.nis = nis_rep;
% ---- 5. combined estimate (common frame)
x = zeros(6,1); XCp = cell(M,1);
for j = 1:M, XCp{j} = local_to_common(S.imm(j).type, S.imm(j).x, S.imm(j).P); x = x + mu(j)*XCp{j}; end
Pc_comb = zeros(6);
for j = 1:M
    [~, Pcj] = local_to_common(S.imm(j).type, S.imm(j).x, S.imm(j).P);
    d = XCp{j} - x; Pc_comb = Pc_comb + mu(j)*(Pcj + d*d');
end
S.x = x; S.P = (Pc_comb + Pc_comb')/2;
end

% ---- native model predict (CV/CA are 6-state; CT is 5-state) ----------------
function [x, Pc] = local_model_predict_n(ty, x, Pc, T, P)
switch ty
    case 'cv'
        F = local_Fcv(T); Q = local_Qcv(T, P.imm_q_cv);
        x = F*x; Pc = F*Pc*F' + Q;
    case 'ca'
        F = local_F(T);   Q = local_Q(T, P.imm_q_ca);
        x = F*x; Pc = F*Pc*F' + Q;
    case 'ct'
        F = local_Jnum(@(z) local_fctn(z, T), x);   % numerical Jacobian of the CT flow
        Q = local_Qct(T, P.imm_q_ct, P.imm_q_omega);
        x = local_fctn(x, T); Pc = F*Pc*F' + Q;
end
end

% ---- native measurement update + Gaussian likelihood ------------------------
function [x, Pc, nis, Lk] = local_update_n(ty, x, Pc, r, th, P)
[R, z] = local_meas(r, th, P);
if strcmp(ty,'ct'), H = [eye(2), zeros(2,3)]; else, H = [eye(2), zeros(2,4)]; end
nu = z - H*x;  Sm = H*Pc*H' + R;  K = (Pc*H')/Sm;
x  = x + K*nu;  IKH = eye(size(Pc,1)) - K*H;
Pc = IKH*Pc*IKH' + K*R*K';  Pc = (Pc + Pc')/2;
nis = (nu'/Sm)*nu;
Lk  = exp(-0.5*nis)/sqrt(max(det(2*pi*Sm), 1e-300));
end

% ---- native prediction over a horizon (no covariance) -----------------------
function xn = local_model_fn(ty, x, h)
switch ty
    case 'cv', xn = local_Fcv(h)*x;
    case 'ca', xn = local_F(h)*x;
    case 'ct', xn = local_fctn(x, h);
end
end

% ---- CT native flow: state [x y vx vy omega], exact coordinated turn --------
function xo = local_fctn(x, T)
px=x(1); py=x(2); vx=x(3); vy=x(4); w=x(5);
if abs(w) > 1e-5
    s = sin(w*T); c = cos(w*T);
    pxn = px + (vx*s - vy*(1-c))/w;
    pyn = py + (vx*(1-c) + vy*s)/w;
    vxn = vx*c - vy*s;
    vyn = vx*s + vy*c;
else
    pxn = px + vx*T; pyn = py + vy*T; vxn = vx; vyn = vy;
end
xo = [pxn; pyn; vxn; vyn; w];
end

function Q = local_Qct(T, qv, qw)
% white-acceleration process noise on (x,vx) and (y,vy) + turn-rate random walk
q2 = qv*[T^3/3 T^2/2; T^2/2 T];
Q = zeros(5);
Q([1 3],[1 3]) = q2;  Q([2 4],[2 4]) = q2;  Q(5,5) = qw*T;
end

% ---- common<->native state/covariance maps (CV,CA identity; CT 5<->6) --------
function [xc, Pc] = local_to_common(ty, x, P)
if strcmp(ty,'ct')
    xc = local_ct2c(x);  J = local_Jnum(@local_ct2c, x);
    if nargout > 1, Pc = J*P*J'; end
else
    xc = x;  if nargout > 1, Pc = P; end
end
end

function [xn, Pn] = local_from_common(ty, xc, Pc)
if strcmp(ty,'ct')
    xn = local_c2ct(xc);  J = local_Jnum(@local_c2ct, xc);
    if nargout > 1, Pn = J*Pc*J'; end
else
    xn = xc;  if nargout > 1, Pn = Pc; end
end
end

function c = local_ct2c(x)                           % CT 5-state -> common 6-state
w = x(5);
c = [x(1); x(2); x(3); x(4); -w*x(4); w*x(3)];       % ax=-w*vy, ay=w*vx (centripetal)
end

function x = local_c2ct(c)                            % common 6-state -> CT 5-state
vx = c(3); vy = c(4); sp = vx*vx + vy*vy;
w = (vx*c(6) - vy*c(5)) / max(sp, 1e-6);             % omega = (vx*ay - vy*ax)/|v|^2
x = [c(1); c(2); vx; vy; w];
end

function J = local_Jnum(f, x)
% Central-difference Jacobian of f: R^n -> R^m (m inferred).
n = numel(x); f0 = f(x); m = numel(f0); J = zeros(m, n);
for k = 1:n
    hk = 1e-6*max(1, abs(x(k)));
    xp = x; xp(k) = xp(k) + hk;
    xm = x; xm(k) = xm(k) - hk;
    J(:,k) = (f(xp) - f(xm))/(2*hk);
end
end

% =========================================================================
% Reference kinematics:  theta_ref and its derivatives from a 6-state vector
% =========================================================================
function [th, thd, thdd] = local_refkin(s)
x = s(1); y = s(2); vx = s(3); vy = s(4); ax = s(5); ay = s(6);
r2 = x*x + y*y; if r2 < 1e-9, r2 = 1e-9; end
th  = atan2(y, x);
n   = x*vy - y*vx;
thd = n / r2;
ndot  = x*ay - y*ax;
r2dot = 2*(x*vx + y*vy);
thdd = ndot/r2 - n*r2dot/(r2*r2);
end

% =========================================================================
% Linear-Gaussian building blocks (CA baseline path is byte-for-byte legacy)
% =========================================================================
function F = local_F(T)                           % constant-acceleration transition
F1 = [1 T T^2/2; 0 1 T; 0 0 1];
F  = kron(F1, eye(2));                            % state order [x y vx vy ax ay]
end

function Q = local_Q(T, q)                        % white-jerk process noise
Q1 = [T^5/20 T^4/8 T^3/6; T^4/8 T^3/3 T^2/2; T^3/6 T^2/2 T];
Q  = q * kron(Q1, eye(2));
end

function F = local_Fcv(T)                          % constant-velocity transition (accel pinned to 0)
F1 = [1 T 0; 0 1 0; 0 0 0];
F  = kron(F1, eye(2));
end

function Q = local_Qcv(T, q)                       % white-acceleration process noise (CV)
Q1 = [T^3/3 T^2/2 0; T^2/2 T 0; 0 0 1e-3];         % small accel-state floor keeps P positive-definite
Q  = q * kron(Q1, eye(2));
end

function [R, z] = local_meas(r, th, P)           % polar -> Cartesian + covariance
c = cos(th); s = sin(th);
z = [r*c; r*s];
J = [c, -r*s; s, r*c];
R = J * diag([P.sig_r^2, P.sig_th^2]) * J';
end

function [x, Pc] = local_kf_init(r, th, P)
[R, z] = local_meas(r, th, P);
x  = [z; 0; 0; 0; 0];
Pc = blkdiag(R, P.sig_v0^2*eye(2), P.sig_a0^2*eye(2));
end

function [x, Pc] = local_kf_predict(x, Pc, dt, P, q)
if nargin < 5, q = P.q_jerk; end
F  = local_F(dt);
x  = F*x;
Pc = F*Pc*F' + local_Q(dt, q);
end

function [x, Pc, nis, Lk] = local_kf_update(x, Pc, r, th, P)
[R, z] = local_meas(r, th, P);
H  = [eye(2), zeros(2,4)];
nu = z - H*x;
Sm = H*Pc*H' + R;
K  = (Pc*H') / Sm;
x  = x + K*nu;
IKH = eye(6) - K*H;
Pc = IKH*Pc*IKH' + K*R*K';                       % Joseph form
Pc = (Pc + Pc')/2;
if nargout > 2
    nis = (nu' / Sm) * nu;                        % normalized innovation squared (2 dof)
end
if nargout > 3
    Lk = exp(-0.5*nis) / sqrt(max(det(2*pi*Sm), 1e-300));   % Gaussian measurement likelihood
end
end
