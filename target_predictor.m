function pred = target_predictor(est, Tp)
%TARGET_PREDICTOR  Predict the target state Tp seconds ahead of the estimate.
%
%   pred = target_predictor(est)       uses Tp = params().t_p (0.5 s)
%   pred = target_predictor(est, Tp)
%
%   Pure function of the estimator struct: no hidden state, no ground truth.
%   Extrapolates the estimated state with the same constant-acceleration model
%   used in the estimator, from est.t to est.t + Tp, and converts the predicted
%   position to elevation as seen from the turret/radar origin.
%
%   pred fields: valid, t_now, Tp, t_target, s (6x1), P (6x6), x_p, y_p, vx_p,
%   vy_p, ax_p, ay_p, r_p, theta_p [rad], sigma_theta [rad, 1-sigma].
%   If est is not initialised, valid = false and all values are NaN.

prm = params();
if nargin < 2 || isempty(Tp)
    Tp = prm.t_p;
end

pred.valid       = false;
pred.t_now       = NaN;
pred.Tp          = Tp;
pred.t_target    = NaN;
pred.s           = nan(6,1);
pred.P           = nan(6,6);
pred.x_p = NaN;  pred.y_p = NaN;  pred.vx_p = NaN;  pred.vy_p = NaN;
pred.ax_p = NaN; pred.ay_p = NaN;
pred.r_p = NaN;  pred.theta_p = NaN;  pred.sigma_theta = NaN;

if ~est.initialized || any(~isfinite(est.s)) || any(~isfinite(est.P(:)))
    return;
end

[F, Q] = ca_model(Tp, est.q_jerk);
s_p = F * est.s;
P_p = F * est.P * F' + Q;
P_p = 0.5 * (P_p + P_p');

pred.t_now    = est.t;
pred.t_target = est.t + Tp;
pred.s        = s_p;
pred.P        = P_p;
pred.x_p  = s_p(1);  pred.y_p  = s_p(2);
pred.vx_p = s_p(3);  pred.vy_p = s_p(4);
pred.ax_p = s_p(5);  pred.ay_p = s_p(6);

pred.r_p     = hypot(pred.x_p, pred.y_p);
pred.theta_p = atan2(pred.y_p, pred.x_p);       % predicted elevation from origin
if pred.r_p > 0
    g = [-pred.y_p, pred.x_p] / pred.r_p^2;     % d(theta)/d(x,y)
    pred.sigma_theta = sqrt(g * P_p(1:2,1:2) * g');
end
pred.valid = true;
end

function [F, Q] = ca_model(dt, q)
% MUST match ca_model in target_state_estimator.m (checked in the test script).
I2 = eye(2);  Z2 = zeros(2);
F = [I2, dt*I2, 0.5*dt^2*I2;
    Z2, I2,    dt*I2;
    Z2, Z2,    I2];
Q1 = q * [dt^5/20, dt^4/8, dt^3/6;
    dt^4/8,  dt^3/3, dt^2/2;
    dt^3/6,  dt^2/2, dt];
Q = zeros(6);
Q([1 3 5],[1 3 5]) = Q1;
Q([2 4 6],[2 4 6]) = Q1;
end