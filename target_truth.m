function [x, y, vx, vy, ax, ay, sid] = target_truth(t, scenario)
%TARGET_TRUTH  Ground-truth target kinematics for scenarios A, B, C.
%   FOR TESTING AND EVALUATION ONLY. The controller must never call this.
%   t may be a scalar or a vector [s]. scenario: 'A','B','C' (any case) or 1,2,3.
%   Outputs: position [m], velocity [m/s], acceleration [m/s^2], scenario id.
%   Frame: turret/radar at origin, x horizontal, y up.

P   = params();
sid = parse_scenario(scenario);
S   = scenario_definition();

% Base motion: constant velocity (Scenario A)
x  = S.x0 + S.vx0 * t;
y  = S.y0 + S.vy0 * t;
vx = S.vx0 * ones(size(t));
vy = S.vy0 * ones(size(t));
ax = zeros(size(t));
ay = zeros(size(t));

% Maneuver (Scenario B, and Scenario C if S.C_has_maneuver)
if sid == 2 || (sid == 3 && S.C_has_maneuver)
    tau = t - P.scnB_maneuver_start;      % time since maneuver onset

    % x: one acceleration bump
    [a1, v1, p1] = unit_bump(tau, S.Tx);
    x  = x  + S.Ax * p1;
    vx = vx + S.Ax * v1;
    ax = ax + S.Ax * a1;

    % y: lane change = positive bump followed by negative bump
    [a2, v2, p2] = unit_bump(tau,          S.Ty);
    [a3, v3, p3] = unit_bump(tau - S.Ty,   S.Ty);
    y  = y  + S.Ay * (p2 - p3);
    vy = vy + S.Ay * (v2 - v3);
    ay = ay + S.Ay * (a2 - a3);
end
end

function S = scenario_definition()
% Engineering assumptions (not given in the assignment PDF)
S.x0 = 1800;  S.y0 = 600;          % initial position [m]
S.vx0 = -100; S.vy0 = 0;           % initial velocity [m/s]
S.Ax = 10;    S.Tx = 4.0;          % x bump: amplitude [m/s^2], duration [s]
S.Ay = 14;    S.Ty = 2.5;          % y half-bump amplitude, duration of each half
S.C_has_maneuver = true;           % Scenario C uses the Scenario B trajectory
end

function [a, v, p] = unit_bump(tau, T)
% Unit-amplitude raised-cosine acceleration pulse and its integrals.
w = 2*pi / T;
a = zeros(size(tau));  v = a;  p = a;

in = (tau >= 0) & (tau <= T);
ti = tau(in);
a(in) = (1 - cos(w*ti)) / 2;
v(in) = (ti - sin(w*ti)/w) / 2;
p(in) = (ti.^2/2 + (cos(w*ti) - 1)/w^2) / 2;

post = tau > T;
tq = tau(post);
v(post) = T/2;
p(post) = T^2/4 + (T/2) * (tq - T);
end

function sid = parse_scenario(scenario)
if ischar(scenario) || (isstring(scenario) && isscalar(scenario))
    s = upper(strtrim(char(scenario)));
    switch s
        case {'A','1'}, sid = 1;
        case {'B','2'}, sid = 2;
        case {'C','3'}, sid = 3;
        otherwise
            error('target_truth:badScenario', 'Unknown scenario "%s" (use A, B or C).', s);
    end
elseif isnumeric(scenario) && isscalar(scenario) && any(scenario == [1 2 3])
    sid = scenario;
else
    error('target_truth:badScenario', 'Scenario must be A/B/C or 1/2/3.');
end
end