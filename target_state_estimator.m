function est = target_state_estimator(est, t, r_meas, theta_meas, meas_valid)
%TARGET_STATE_ESTIMATOR  Extended Kalman filter for the target (elevation plane).
%
%   est = target_state_estimator()
%       Fresh, uninitialised estimator (also how to reset).
%   est = target_state_estimator(est, t, r_meas, theta_meas, meas_valid)
%       One step: propagate to time t, then update if a NEW valid radar
%       measurement is supplied. Call at every controller step.
%
%   State s = [x y vx vy ax ay]'  [m, m, m/s, m/s, m/s^2, m/s^2], origin at the
%   turret/radar, x horizontal, y up. Motion model: constant acceleration with
%   white-noise jerk (spectral density est.q_jerk). Measurement: range and
%   elevation, h(s) = [hypot(x,y); atan2(y,x)].
%
%   Uses only radar measurements and its own past estimate. It never reads
%   ground truth and has no knowledge of scenarios or blackout timing.
%
%   Fields of est:
%     initialized, s (6x1), P (6x6), t (time of s), t_init, t_meas (time of last
%     measurement used), n_updates (measurements used, incl. the initialising
%     one), updated (true if a measurement was used in this call), innov (2x1),
%     nis (normalised innovation squared), and tuning: sigma_r, sigma_theta,
%     sigma_v0, sigma_a0, q_jerk.

    if nargin == 0
        est = new_estimator();
        return;
    end
    if isempty(est)
        est = new_estimator();
    end

    % Time moved backwards => a new run has started: restart, keep tuning.
    if est.initialized && (t < est.t - 1e-9)
        est = restart_keep_tuning(est);
    end

    have_meas   = (meas_valid ~= 0) && isfinite(r_meas) && isfinite(theta_meas);
    est.updated = false;

    % ---- Not initialised yet: wait for the first valid measurement ----
    if ~est.initialized
        if have_meas
            est = initialise(est, t, r_meas, theta_meas);
        end
        return;
    end

    % ---- Time update (prediction step only) ----
    dt = t - est.t;
    if dt > 0
        [F, Q] = ca_model(dt, est.q_jerk);
        est.s = F * est.s;
        est.P = F * est.P * F' + Q;
        est.P = 0.5 * (est.P + est.P');
        est.t = t;
    end

    % ---- Measurement update (only for a new valid measurement) ----
    if have_meas && (t - est.t_meas > 1e-9)
        est = ekf_update(est, r_meas, theta_meas);
        if est.updated
            est.t_meas = t;
        end
    end
end

% =====================================================================
function est = new_estimator()
    prm = params();
    est.initialized = false;
    est.s           = nan(6,1);
    est.P           = nan(6,6);
    est.t           = NaN;
    est.t_init      = NaN;
    est.t_meas      = NaN;
    est.n_updates   = 0;
    est.updated     = false;
    est.innov       = nan(2,1);
    est.nis         = NaN;

    % ---- Tuning / design parameters ----
    est.sigma_r     = prm.sigma_r;          % range noise std [m]
    est.sigma_theta = prm.sigma_theta_rad;  % elevation noise std [rad]
    est.sigma_v0    = prm.v_max / 2;        % prior velocity std per axis [m/s]
    est.sigma_a0    = prm.a_max / 2;        % prior acceleration std per axis [m/s^2]
    est.q_jerk      = 10;                   % white-jerk PSD [m^2/s^5] (see design notes)
end

function est = restart_keep_tuning(old)
    est = new_estimator();
    names = {'sigma_r','sigma_theta','sigma_v0','sigma_a0','q_jerk'};
    for k = 1:numel(names)
        est.(names{k}) = old.(names{k});
    end
end

function est = initialise(est, t, r_m, th_m)
    c  = cos(th_m);
    sn = sin(th_m);
    est.s = [r_m*c; r_m*sn; 0; 0; 0; 0];

    J  = [c, -r_m*sn; sn, r_m*c];                      % d(x,y)/d(r,theta)
    Rp = diag([est.sigma_r^2, est.sigma_theta^2]);
    est.P = zeros(6);
    est.P(1:2,1:2) = J * Rp * J';
    est.P(3:4,3:4) = est.sigma_v0^2 * eye(2);
    est.P(5:6,5:6) = est.sigma_a0^2 * eye(2);

    est.initialized = true;
    est.t           = t;
    est.t_init      = t;
    est.t_meas      = t;
    est.n_updates   = 1;
    est.updated     = true;
end

function est = ekf_update(est, r_m, th_m)
    x = est.s(1);  y = est.s(2);
    r_pred = hypot(x, y);
    if r_pred < 1e-6
        return;                       % degenerate geometry: skip update
    end
    th_pred = atan2(y, x);

    H = zeros(2,6);
    H(1,1) =  x / r_pred;      H(1,2) = y / r_pred;
    H(2,1) = -y / r_pred^2;    H(2,2) = x / r_pred^2;
    R = diag([est.sigma_r^2, est.sigma_theta^2]);

    nu = [r_m - r_pred; wrap_angle(th_m - th_pred)];
    S  = H * est.P * H' + R;
    K  = (est.P * H') / S;

    est.s = est.s + K * nu;
    IKH   = eye(6) - K * H;
    est.P = IKH * est.P * IKH' + K * R * K';      % Joseph form
    est.P = 0.5 * (est.P + est.P');

    est.innov     = nu;
    est.nis       = nu' * (S \ nu);
    est.n_updates = est.n_updates + 1;
    est.updated   = true;
end

function [F, Q] = ca_model(dt, q)
    % Constant-acceleration model, white-noise jerk. State [x y vx vy ax ay].
    % (Same model is used in target_predictor.m.)
    I2 = eye(2);  Z2 = zeros(2);
    F = [I2, dt*I2, 0.5*dt^2*I2;
         Z2, I2,    dt*I2;
         Z2, Z2,    I2];
    Q1 = q * [dt^5/20, dt^4/8, dt^3/6;
              dt^4/8,  dt^3/3, dt^2/2;
              dt^3/6,  dt^2/2, dt];
    Q = zeros(6);
    Q([1 3 5],[1 3 5]) = Q1;           % x-axis: x, vx, ax
    Q([2 4 6],[2 4 6]) = Q1;           % y-axis: y, vy, ay
end

function a = wrap_angle(a)
    a = atan2(sin(a), cos(a));         % wrap to (-pi, pi]
end