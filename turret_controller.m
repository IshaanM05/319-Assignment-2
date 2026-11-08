function [V_m, fire_cmd, dbg] = turret_controller(t, r_meas, theta_meas, meas_valid, theta_turret, theta_dot_turret) %#ok<INUSD>
%TURRET_CONTROLLER  Predictive tracking controller for the elevation mount.
%
%   Data flow (Sec. 7):  polar measurement -> CA Kalman filter (Cartesian)
%     -> tp = 0.5 s prediction -> theta_ref = atan2(y_pred, x_pred) (clamped)
%     -> LTI compensator Gc(s) = Kp/(1 + s/p)  (+ voltage saturation, anti-windup)
%     -> V_m, plus single-sample fire_cmd logic with a 1 s refractory period.
%
%   Persistent state (cleared with `clear turret_controller`; also re-initialised
%   automatically if time ever decreases): KF state/covariance, compensator
%   state, last pulse time, previous fire_cmd.
%
%   The controller may be called at ANY rate (all discretisations use the
%   actual elapsed time t - t_prev).  theta_dot_turret is part of the
%   required interface but is not used by the final (single-loop) design.
%
%   Optional 3rd output `dbg` (only built if requested) exposes internal
%   signals for plotting.
persistent S
P = turret_params();
if isempty(S) || t < S.t_prev - 1e-9
    S = local_new_state(t);
end
dt = t - S.t_prev;

%% ---------------------------------------------------------------- estimator
if ~S.init
    if meas_valid && isfinite(r_meas) && isfinite(theta_meas)
        [S.x, S.P] = local_kf_init(r_meas, theta_meas, P);
        S.init = true;
    end
else
    if dt > 0, [S.x, S.P] = local_kf_predict(S.x, S.P, dt, P); end
    if meas_valid && isfinite(r_meas) && isfinite(theta_meas)
        [S.x, S.P] = local_kf_update(S.x, S.P, r_meas, theta_meas, P);   % else: propagate only
    end
end
S.t_prev = t;

%% --------------------------------------------- prediction -> elevation reference
if S.init
    if strcmp(P.pred_model, 'CV')
        xp = S.x(1:2) + S.x(3:4)*P.tp;
    else
        xp = local_F(P.tp) * S.x;  xp = xp(1:2);                      % x_pred = F(tp) x
    end
    th_ref = atan2(xp(2), xp(1));
    th_ref = min(max(th_ref, P.th_min), P.th_max);                    % mount limits
else
    xp = [NaN; NaN];  th_ref = theta_turret;                          % nothing to track yet
end

%% ----------------------------------------------- LTI compensator + saturation
e = th_ref - theta_turret;
if ~S.init
    V_unsat = 0; V_m = 0;
elseif strcmp(P.ctrl, 'pi')                       % ablation: suggested PI
    xi_new = S.xi + e*dt;
    V_unsat = P.Kp_pi*(e + xi_new/P.Ti_pi);
    V_m = min(max(V_unsat, -P.Vmax), P.Vmax);
    if V_m ~= V_unsat                              % back-calculation anti-windup
        xi_new = P.Ti_pi*(V_m/P.Kp_pi - e);
    end
    S.xi = xi_new;
else                                              % FINAL: Gc(s) = Kp/(1+s/p)
    a = 1 - exp(-P.p_eff*dt);                      % exact first-order discretisation
    xf = S.xf + a*(P.Kp_eff*e - S.xf);
    V_unsat = xf;
    V_m = min(max(xf, -P.Vmax), P.Vmax);
    S.xf = V_m;                                    % state == output -> clamping state = exact anti-windup
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
    dbg.x = S.x; dbg.sigma = sqrt(max(diag(S.P), 0)); dbg.xp = xp;
    dbg.th_ref = th_ref; dbg.V_unsat = V_unsat; dbg.dperp = dperp; dbg.sig_perp = sig_perp;
    dbg.init = S.init;
end
end

% =========================================================================
function S = local_new_state(t)
S.init = false; S.x = zeros(6,1); S.P = eye(6);
S.t_prev = t; S.xf = 0; S.xi = 0;
S.t_fire = -inf; S.fire_prev = false;
end

function F = local_F(T)                           % constant-acceleration transition
F1 = [1 T T^2/2; 0 1 T; 0 0 1];
F  = kron(F1, eye(2));                            % state order [x y vx vy ax ay]
end

function Q = local_Q(T, q)                        % white-jerk process noise
Q1 = [T^5/20 T^4/8 T^3/6; T^4/8 T^3/3 T^2/2; T^3/6 T^2/2 T];
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

function [x, Pc] = local_kf_predict(x, Pc, dt, P)
F  = local_F(dt);
x  = F*x;
Pc = F*Pc*F' + local_Q(dt, P.q_jerk);
end

function [x, Pc] = local_kf_update(x, Pc, r, th, P)
[R, z] = local_meas(r, th, P);
H  = [eye(2), zeros(2,4)];
Sm = H*Pc*H' + R;
K  = (Pc*H') / Sm;
x  = x + K*(z - H*x);
IKH = eye(6) - K*H;
Pc = IKH*Pc*IKH' + K*R*K';                       % Joseph form
Pc = (Pc + Pc')/2;
end
