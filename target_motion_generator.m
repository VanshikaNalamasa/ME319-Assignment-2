function [r_meas, theta_meas, meas_valid] = target_motion_generator(t, scenario)
%TARGET_MOTION_GENERATOR  Radar measurements of the V-44 target.
%   t        : time [s] (scalar)
%   scenario : 'A','B','C' (any case) or 1,2,3
%   r_meas, theta_meas : noisy range [m] and elevation [rad]; NaN when no
%                        new measurement is available
%   meas_valid         : true only at a 20 Hz radar instant outside a blackout
%
%   Stateless: the output depends only on (t, scenario). Noise is seeded by
%   the radar sample index, so it is reproducible, identical on repeated calls
%   and independent of call order. Global random state is not touched.
%   True target state is never returned.

P = params();

t_tol     = 1e-6;        % [s] tolerance for "t is on a radar instant"
seed_base = 20440000;    % fixed seed offset for reproducible noise

% Default: no measurement
r_meas     = NaN;
theta_meas = NaN;
meas_valid = false;

% Nearest radar sample instant
k   = round(t / P.Ts_radar);
t_s = k * P.Ts_radar;

% Target position at that instant (also validates 'scenario')
[x, y, ~, ~, ~, ~, sid] = target_truth(t_s, scenario);

% New measurement only exactly on a 20 Hz instant
if ~((k >= 0) && (abs(t - t_s) <= t_tol))
    return;
end

% Scenario C blackout window [t_bo0, t_bo0 + t_bo)
t0 = P.scnC_t_bo0;
t1 = P.scnC_t_bo0 + P.scnC_t_bo;
if sid == 3 && (t_s >= t0 - t_tol) && (t_s < t1 - t_tol)
    return;
end

% Noisy polar measurement (independent Gaussian noise)
stream = RandStream('mt19937ar', 'Seed', seed_base + k);
n = randn(stream, 2, 1);
r_meas     = hypot(x, y) + P.sigma_r         * n(1);
theta_meas = atan2(y, x) + P.sigma_theta_rad * n(2);
meas_valid = true;
end