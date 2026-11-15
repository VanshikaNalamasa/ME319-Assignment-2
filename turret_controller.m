function [V_m, fire_cmd] = turret_controller(t, r_meas, theta_meas, meas_valid, theta_turret, theta_dot_turret) %#ok<INUSD>
%TURRET_CONTROLLER  Predictive-feedback controller for the turret elevation axis.
%
%   [V_m, fire_cmd] = turret_controller(t, r_meas, theta_meas, meas_valid, ...
%                                       theta_turret, theta_dot_turret)
%
%   Call ONCE PER SIMULATION STEP with the radar outputs of that step.
%
%   ARCHITECTURE
%     radar -> target_state_estimator (EKF) -> target_predictor (Tp = params().t_p)
%           -> theta_ref = clamp(predicted elevation, [-10 deg, 85 deg])
%           -> V_m = Kp*(theta_ref - theta_turret)          (LTI proportional loop)
%     Laser firing logic is SEPARATE (fire_logic below) and never feeds back
%     into V_m.
%
%   FEEDBACK DESIGN (all from the plant in params.m)
%     Gp(s) = G_num / ( s*(G_den(1)s^2 + G_den(2)s + G_den(3)) )
%     Low-frequency velocity constant  Kv = G_num/G_den(3)   [rad/s per V]
%     Gain                             Kp = 1/(Tp*Kv)
%     => loop velocity constant Kp*Kv = 1/Tp (type-1 loop).
%     For a RAMP reference the steady-state lag is 1/(Kp*Kv) = Tp seconds,
%     so the turret lags the predicted reference (theta(t+Tp)) by about Tp,
%     i.e. it points near theta(t).  This holds ONLY in steady state for a
%     constant-rate reference.  Transients (initial slew, maneuver onset,
%     measurement recovery) and angular acceleration are NOT cancelled:
%     the quasi-static residual pointing error is ~ c*theta_ddot with
%     c = Tp^2/2 - c1*Tp + (c1^2 - c2), c1 = G_den(3)/(G_num*Kp),
%     c2 = G_den(2)/(G_num*Kp)   (derived in the test script and checked there).
%
%   FIRE LOGIC (separate from feedback)
%     The beam is a ray at the instantaneous turret angle and a hit is judged
%     against the target position at the firing instant.  So the gate uses the
%     ESTIMATED CURRENT target position p_hat (not the prediction):
%        n          = [-sin(theta_turret); cos(theta_turret)]   (normal to beam)
%        d_hat      = |n' * p_hat|                              (miss distance)
%        sigma_perp = sqrt( n' * P_xy * n )                     (1-sigma, same direction)
%        gate open  <=>  target in front of beam  AND  d_hat + k_sigma*sigma_perp <= R_eng
%     Fire only if the gate is open AND at least t_recharge (+1e-6 s guard against
%     floating-point round-off) has elapsed since the previous pulse.
%     fire_cmd is HIGH FOR EXACTLY ONE STEP per pulse, then forced LOW, so every
%     pulse is a clean 0->1 transition.
%
%   ASSUMPTIONS (state them in the report)
%     A1  Controller is called at every step of a fixed-step loop, t = k*dt.
%     A2  The laser pulse lasts one simulation step (turret angle constant during
%         it); the pulse is evaluated at the step where fire_cmd goes 0->1.
%     A3  No voltage/current saturation is modelled (none given in the assignment).
%     A4  State (estimator, recharge timer) is kept in a persistent variable
%         because the required interface has no state argument.  It is reset
%         whenever t decreases (a new run) or when turret_controller() is called
%         with no arguments.
%     A5  Turret limits are respected through the REFERENCE clamp.
%
%   SERVICE CALLS (used by the test script only)
%     turret_controller()          reset internal state
%     turret_controller('reset')   same
%     turret_controller('design')  returns struct with Kv, Kp, Tp, Kv_loop, k_sigma, ...
%     turret_controller('info')    returns diagnostics of the most recent step
%
%   Uses ONLY radar inputs and the turret state. No ground truth.

persistent st

% ------------------------------------------------------------ service calls
if nargin == 0
    st = [];  V_m = 0;  fire_cmd = false;
    return;
end
if nargin == 1 && ischar(t)
    switch lower(t)
        case 'reset'
            st = [];  V_m = 0;
        case 'design'
            V_m = controller_design();
        case 'info'
            if isempty(st), V_m = struct(); else, V_m = st.info; end
        otherwise
            error('turret_controller:badCommand', 'Unknown command "%s".', t);
    end
    fire_cmd = false;
    return;
end

% ------------------------------------------------- (re)initialise on new run
if isempty(st) || t < st.last_t - 1e-9
    st = new_state();
end
st.last_t = t;
prm = st.prm;
des = st.des;

% ------------------------------------------------------- estimate + predict
st.est = target_state_estimator(st.est, t, r_meas, theta_meas, meas_valid);
pred   = target_predictor(st.est, prm.t_p);

info = blank_info(t);
info.est_n_updates = st.est.n_updates;
info.est_t_init    = st.est.t_init;
info.est_updated   = st.est.updated;

if ~pred.valid                      % estimator not initialised yet
    V_m = 0;  fire_cmd = false;  st.prev_fire = false;
    st.info = info;
    return;
end

% -------------------------------------------------------- feedback (LTI, P)
theta_ref = min(max(pred.theta_p, prm.theta_min_rad), prm.theta_max_rad);
V_m = des.Kp * (theta_ref - theta_turret);
if ~isfinite(V_m), V_m = 0; end

% -------------------------------------------------- firing logic (separate)
[fire_cmd, g] = fire_logic(st, t, theta_turret);
st.prev_fire = fire_cmd;
if fire_cmd
    st.last_pulse_t = t;
end

% ------------------------------------------------------------- diagnostics
info.valid         = true;
info.theta_ref     = theta_ref;
info.theta_ref_raw = pred.theta_p;
info.V_m           = V_m;
info.d_hat         = g.d_hat;
info.sigma_perp    = g.sigma_perp;
info.along         = g.along;
info.gate_open     = g.gate_open;
info.recharge_ok   = g.recharge_ok;
info.fire          = fire_cmd;
info.theta_hat     = atan2(st.est.s(2), st.est.s(1));
st.info = info;
end

% =====================================================================
function des = controller_design()
    prm = params();
    des.Kv             = prm.G_num / prm.G_den(3);   % rad/s per V (low-frequency)
    des.Tp             = prm.t_p;
    des.Kp             = 1 / (des.Tp * des.Kv);      % V/rad
    des.Kv_loop        = des.Kp * des.Kv;            % = 1/Tp
    des.k_sigma        = 1.0;                        % gate: d_hat + k*sigma_perp <= R_eng
    des.recharge_guard = 1e-6;                       % [s] round-off guard on 1.0 s recharge
end

function st = new_state()
    st.prm          = params();
    st.des          = controller_design();
    st.est          = target_state_estimator();
    st.last_t       = -Inf;
    st.last_pulse_t = -Inf;
    st.prev_fire    = false;
    st.info         = blank_info(NaN);
end

function info = blank_info(t)
    info.t             = t;
    info.valid         = false;
    info.theta_ref     = NaN;
    info.theta_ref_raw = NaN;
    info.V_m           = 0;
    info.d_hat         = NaN;
    info.sigma_perp    = NaN;
    info.along         = NaN;
    info.gate_open     = false;
    info.recharge_ok   = false;
    info.fire          = false;
    info.theta_hat     = NaN;
    info.est_n_updates = 0;
    info.est_t_init    = NaN;
    info.est_updated   = false;
end

function [fire, g] = fire_logic(st, t, theta_t)
    est = st.est;  prm = st.prm;  des = st.des;

    p = est.s(1:2);                          % estimated CURRENT target position
    u = [cos(theta_t);  sin(theta_t)];       % beam direction
    n = [-sin(theta_t); cos(theta_t)];       % unit normal to the beam
    along = u' * p;                          % distance along the beam

    Pxy        = est.P(1:2, 1:2);
    sigma_perp = sqrt(max(n' * Pxy * n, 0)); % 1-sigma position uncertainty normal to beam

    if along > 0
        d_hat = abs(n' * p);                 % distance from target to the beam ray
    else
        d_hat = norm(p);                     % target behind turret: closest point is origin
    end

    gate_open   = (along > 0) && (d_hat + des.k_sigma * sigma_perp <= prm.R_eng);
    recharge_ok = (t - st.last_pulse_t) >= (prm.t_recharge + des.recharge_guard);

    fire = logical(gate_open && recharge_ok && ~st.prev_fire);

    g.d_hat = d_hat;  g.sigma_perp = sigma_perp;  g.along = along;
    g.gate_open = gate_open;  g.recharge_ok = recharge_ok;
end
