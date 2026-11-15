function [theta_dot_turret, theta_ddot_turret, i_dot] = ...
    turret_plant(t, V_m, theta_turret, theta_dot_turret, i) %#ok<INUSL>
%TURRET_PLANT  State derivatives of the elevation drive (motor + gear + turret).
%   Inputs : t [s] (unused, time-invariant), V_m [V], theta_turret [rad],
%            theta_dot_turret [rad/s], i [A]
%   Outputs: theta_dot_turret [rad/s] (returned unchanged),
%            theta_ddot_turret [rad/s^2], i_dot [A/s]
%
%   Self-contained: parameters are hard-coded (same values as params.m).

% Parameters
J_turret = 1200;    % [kg m^2]
J_motor  = 0.20;    % [kg m^2]
N        = 100;     % gear ratio
b        = 1500;    % [N m s/rad], turret side
Kt       = 10;      % [N m/A]
Ke       = 10;      % [V s/rad]
Rm       = 1.0;     % [ohm]
Lm       = 0.10;    % [H]

J_eq = J_turret + N^2 * J_motor;    % 3200 kg m^2

% Mechanical: J_eq*theta_ddot = N*Kt*i - b*theta_dot
theta_ddot_turret = (N*Kt*i - b*theta_dot_turret) / J_eq;

% Electrical: Lm*di/dt = V - Rm*i - Ke*N*theta_dot
i_dot = (V_m - Rm*i - N*Ke*theta_dot_turret) / Lm;

% theta_dot_turret is returned unchanged (same name as the input).
end