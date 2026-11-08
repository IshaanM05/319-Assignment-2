function [r_meas, theta_meas, meas_valid, truth] = target_motion_generator(t, scenario)
%TARGET_MOTION_GENERATOR  Object trajectory + noisy 20 Hz polar radar measurements.
%
%   [r_meas, theta_meas, meas_valid] = target_motion_generator(t, scenario)
%   [..., truth] = ...   optional 4th output: true state at time t
%                        (x,y,vx,vy,ax,ay,in_loss) - used by the simulation
%                        driver for scoring ONLY, never by the controller.
%
%   scenario : 'A' non-manoeuvring | 'B' manoeuvring (from t = 1 s)
%              'C' = 'B' + sensor loss of 0.50 s starting at t = 1.00 s
%              (numeric 1/2/3 are accepted as aliases)
%              struct with fields  .seed  (integer: randomised bounded manoeuvre,
%              noise realisation derived from the seed) and optionally
%              .loss = [t_onset dur; ...]  (rows; each dur <= 0.5, onset >= 1).
%
%   Measurements are produced when t crosses a multiple of 0.05 s; every other
%   call returns meas_valid = false and NaN measurements.  Samples inside a loss
%   interval [t0, t0+dur) are dropped.
%
%   Envelope enforced for EVERY generated trajectory (error / rejection otherwise):
%       800 <= r <= 2000 m,  |v| <= 150 m/s,  |a| <= 20 m/s^2,  smooth a(t).
%   Clear state between runs with  `clear target_motion_generator`.
persistent C
if nargin < 2 || isempty(scenario), scenario = 'A'; end
key = local_key(scenario);
if isempty(C) || ~strcmp(C.key, key) || t < C.t_last - 1e-9
    C = local_build(scenario);  C.key = key;  C.last_k = -1;
end
C.t_last = t;

r_meas = NaN; theta_meas = NaN; meas_valid = false;
k = floor(t/C.dtm + 1e-9);
if k > C.last_k
    C.last_k = k;
    tk = k*C.dtm;
    if (k+1 <= numel(C.nr)) && ~local_in_loss(C.loss, tk)
        s = local_truth(C, tk);
        r_meas     = hypot(s.x, s.y) + C.nr(k+1);
        theta_meas = atan2(s.y, s.x) + C.nth(k+1);
        meas_valid = true;
    end
end
if nargout > 3
    truth = local_truth(C, t);
    truth.in_loss = local_in_loss(C.loss, t);
end
end

% =========================================================================
function key = local_key(sc)
if isstruct(sc) && isfield(sc,'base')            % named motion + custom loss intervals
    key = [upper(sc.base) mat2str(sc.loss)];
elseif isstruct(sc)
    key = sprintf('S%d', sc.seed);
    if isfield(sc,'loss'), key = [key mat2str(sc.loss)]; end
elseif isnumeric(sc)
    key = char('A' + sc - 1);
else
    key = upper(char(sc));
end
end

function n = local_name(sc)
if isstruct(sc), n = sc.base; elseif isnumeric(sc), n = char('A' + sc - 1); else, n = char(sc); end
end

function tf = local_in_loss(loss, t)
tf = false;
for j = 1:size(loss,1)
    if t >= loss(j,1) - 1e-9 && t < loss(j,1) + loss(j,2) - 1e-9, tf = true; return; end
end
end

function s = local_truth(C, t)
u = min(max(t, 0), C.T)/C.h;  i = min(floor(u), C.n-2);  f = u - i;  i = i + 1;
g = @(v) v(i)*(1-f) + v(i+1)*f;
s.x = g(C.X); s.y = g(C.Y); s.vx = g(C.VX); s.vy = g(C.VY); s.ax = g(C.AX); s.ay = g(C.AY);
end

% -------------------------------------------------------------------------
function C = local_build(sc)
P = turret_params();
C.dtm = P.dt_meas;  C.h = 1e-3;  C.T = 10.5;  C.n = round(C.T/C.h) + 1;
tg = (0:C.n-1)' * C.h;
if isstruct(sc) && ~isfield(sc,'base')
    seed = sc.seed;
    [pos0, v0, segs] = local_random_traj(seed, tg, P);
    if isfield(sc,'loss')
        C.loss = sc.loss;
    else
        C.loss = local_random_loss(seed);
    end
    nseed = 7919*seed + 13;
else
    nm_ = upper(local_name(sc));
    switch nm_
        case 'A'      % non-manoeuvring: constant velocity
            pos0 = [1250; 650];  v0 = [-45; 12];  segs = zeros(0,5);  C.loss = zeros(0,2);
        case {'B','C'}
            pos0 = [1300; 600];  v0 = [-35; -5];
            %        t_start  dur   ax    ay   shape(1 bump, 2 S-turn, 3 hold)
            segs = [ 1.0      5.0  -10.0  15.0   1;
                     5.5      4.0    9.0 -12.0   2 ];
            C.loss = zeros(0,2);
            if strcmp(nm_,'C'), C.loss = [1.00 0.50]; end
            if isstruct(sc), C.loss = sc.loss; end
        otherwise
            error('target_motion_generator:scenario', 'Unknown scenario "%s".', nm_);
    end
    nseed = 12345 + double(nm_);
    if strcmp(nm_,'C'), nseed = 12345 + double('B'); end   % C = same motion AND noise as B
end
[X, Y, VX, VY, AX, AY] = local_integrate(tg, pos0, v0, segs);
ok = local_check(tg, X, Y, VX, VY, AX, AY, P, false);
if ~ok && (~isstruct(sc) || isfield(sc,'base')), error('target_motion_generator:bounds', 'Scenario violates the operating envelope.'); end
C.X = X; C.Y = Y; C.VX = VX; C.VY = VY; C.AX = AX; C.AY = AY;
% measurement-noise realisation (independent of how often the function is called)
rs = RandStream('mt19937ar', 'Seed', nseed);
nm = round(10/C.dtm) + 2;
C.nr  = P.sig_r  * randn(rs, nm, 1);
C.nth = P.sig_th * randn(rs, nm, 1);
end

function [X, Y, VX, VY, AX, AY] = local_integrate(tg, pos0, v0, segs)
AX = zeros(size(tg)); AY = AX;
for j = 1:size(segs,1)
    ts = segs(j,1); D = segs(j,2); sh = segs(j,5);
    u = (tg - ts)/D;
    w = zeros(size(tg));
    switch sh
        case 1, m = (u>=0 & u<=1); w(m) = sin(pi*u(m)).^2;                       % smooth bump
        case 2, m = (u>=0 & u<=1); w(m) = sin(pi*u(m)).^2 .* sin(2*pi*u(m));     % S-turn (zero net dv)
        case 3, q = min(max((tg-ts)/0.5, 0), 1); w = q.^2 .* (3-2*q);            % smooth-step then hold
    end
    AX = AX + segs(j,3)*w;  AY = AY + segs(j,4)*w;
end
h = tg(2) - tg(1);
VX = v0(1) + [0; cumsum((AX(1:end-1)+AX(2:end))/2*h)];
VY = v0(2) + [0; cumsum((AY(1:end-1)+AY(2:end))/2*h)];
X  = pos0(1) + [0; cumsum((VX(1:end-1)+VX(2:end))/2*h)];
Y  = pos0(2) + [0; cumsum((VY(1:end-1)+VY(2:end))/2*h)];
end

function ok = local_check(tg, X, Y, VX, VY, AX, AY, P, strictTheta) %#ok<INUSD>
m = tg <= 10 + 1e-9;
r = hypot(X(m), Y(m)); v = hypot(VX(m), VY(m)); a = hypot(AX(m), AY(m));
th = atan2(Y(m), X(m));
ok = all(r >= P.r_min) && all(r <= P.r_max) && all(v <= P.v_max) && all(a <= P.a_max + 1e-9) ...
     && all(th >= deg2rad(3)) && all(th <= deg2rad(80)) && all(X(m) > 0);
end

function [pos0, v0, segs] = local_random_traj(seed, tg, P)
% Randomised bounded manoeuvre: 1-3 smooth acceleration segments, the first
% starting at t >= 1 s, in a random direction (so both x and y move).
rs = RandStream('mt19937ar', 'Seed', seed);
for attempt = 1:5000
    r0  = P.r_min + 100 + rand(rs)*(P.r_max - P.r_min - 300);   % 900 .. 1900 m
    th0 = deg2rad(5 + rand(rs)*65);                              % 5 .. 70 deg
    pos0 = r0*[cos(th0); sin(th0)];
    sp0 = rand(rs)*100;  hd = 2*pi*rand(rs);
    v0 = sp0*[cos(hd); sin(hd)];
    nseg = randi(rs, 3);
    segs = zeros(nseg, 5);  ts = 1 + 2*rand(rs);
    for j = 1:nseg
        D  = 1.5 + 3*rand(rs);
        A  = 8 + 12*rand(rs);   phi = 2*pi*rand(rs);
        sh = 1 + (rand(rs) > 0.4) + (rand(rs) > 0.7);            % mix of bump / S-turn / hold
        segs(j,:) = [ts, D, A*cos(phi), A*sin(phi), min(sh,3)];
        ts = ts + 0.5 + 2.5*rand(rs);
        if ts > 8, ts = 8; end
    end
    % scale so that the peak |a| respects a_max (with a 5% margin)
    [~,~,~,~,AX,AY] = local_integrate(tg, pos0, v0, segs);
    pk = max(hypot(AX(tg<=10), AY(tg<=10)));
    if pk > 0.95*P.a_max, segs(:,3:4) = segs(:,3:4) * (0.95*P.a_max/pk); end
    [X,Y,VX,VY,AX,AY] = local_integrate(tg, pos0, v0, segs);
    if local_check(tg, X, Y, VX, VY, AX, AY, P, true), return; end   % accept
end
error('target_motion_generator:reject', 'No admissible random trajectory found.');
end

function loss = local_random_loss(seed)
% one or two loss intervals; durations in [0.1, 0.5] s; onset >= 1.0 s
rs = RandStream('mt19937ar', 'Seed', 31*seed + 7);
n = 1 + (rand(rs) < 0.3);
loss = zeros(n,2);
t0 = 1 + 7.5*rand(rs);
for j = 1:n
    d = 0.1 + 0.4*rand(rs);
    loss(j,:) = [t0, d];
    t0 = t0 + d + 0.5 + 1.0*rand(rs);
    if t0 > 9.4, loss = loss(1:j,:); break; end
end
end
