function [theta_dot_turret, theta_ddot_turret, i_dot] = turret_plant(t, V_m, theta_turret, theta_dot_turret, i) %#ok<INUSL>
%TURRET_PLANT  Elevation-drive dynamics of the laser turret (states [theta, theta_dot, i]).
%
%   The elevation drive is a motor coupled to the turret through a 100:1 gear
%   reduction.  Reflecting the rotor inertia through the gearbox gives
%   J_eq = J_turret + N^2*J_m = 3200 kg m^2, and the two first-principles
%   equations of the brief (Sec. 2) are
%
%     mechanical :  J_eq*theta_ddot + b*theta_dot = N*Kt*i
%     electrical :  Lm*di/dt + Rm*i + Ke*N*theta_dot = V_m
%
%   with b the turret-side viscous damping, Rm/Lm the armature resistance and
%   inductance, Kt/Ke the torque and back-emf constants, and N = 100 the gear
%   ratio (the motor turns N times faster than the turret).  Eliminating i(s)
%   gives the reference plant Gp(s) used for the linear design (paper Eq. (3)).
%
%   Returns [d(theta)/dt, d(theta_dot)/dt, di/dt].  The first output equals the
%   input theta_dot_turret except when the turret sits on a mechanical hard stop
%   (-10 deg / +85 deg) and is being driven outward: then the position
%   derivative is zeroed and the velocity is damped out with a 5 ms time
%   constant (stiff mechanical end stop).
persistent Pp
if isempty(Pp), Pp = turret_params(); end
w = theta_dot_turret;                                   % keep the input value
theta_ddot_turret = (Pp.N*Pp.Kt*i - Pp.b*w) / Pp.J_eq;
i_dot             = (V_m - Pp.Rm*i - Pp.Ke*Pp.N*w) / Pp.Lm;
if (theta_turret >= Pp.th_max && w > 0) || (theta_turret <= Pp.th_min && w < 0)
    theta_dot_turret  = 0;                              % position frozen at the stop
    theta_ddot_turret = -w/5e-3;                        % stop absorbs the velocity
end
end
