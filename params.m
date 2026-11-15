function p = params()
%PARAMS  All ME319 Assignment 2 parameters (values from the assignment PDF).
%   p = params() returns a struct. Units are SI unless the field name says otherwise.

%% Turret and motor mechanics (Section 2)
p.J_turret = 1200;      % turret + cradle inertia about trunnion [kg m^2]
p.J_motor  = 0.20;      % motor rotor inertia [kg m^2]
p.N        = 100;       % gear ratio (motor speed = N * turret speed)
p.b        = 1500;      % turret-side viscous damping [N m s/rad]

%% Motor electrical
p.Kt = 10;              % torque constant [N m/A]
p.Ke = 10;              % back-emf constant [V s/rad]
p.Rm = 1.0;             % armature resistance [ohm]
p.Lm = 0.10;            % armature inductance [H]

%% Derived: equivalent inertia reflected to the turret side
p.J_eq = p.J_turret + p.N^2 * p.J_motor;   % = 3200 [kg m^2]

%% Derived: plant Theta(s)/V(s) = G_num / (G_den(1) s^3 + ... + G_den(4) s)
p.G_num = p.N * p.Kt;
p.G_den = [ p.Lm*p.J_eq, ...
    p.Lm*p.b + p.Rm*p.J_eq, ...
    p.Rm*p.b + p.N^2*p.Kt*p.Ke, ...
    0 ];

%% Turret angle limits
p.theta_min_deg = -10;
p.theta_max_deg = 85;
p.theta_min_rad = p.theta_min_deg * pi/180;
p.theta_max_rad = p.theta_max_deg * pi/180;

%% Radar
p.f_radar         = 20;                 % [Hz]
p.Ts_radar        = 1 / p.f_radar;      % [s] (= 0.05)
p.sigma_r         = 1;                  % range noise std [m]
p.sigma_theta_deg = 0.02;               % angle noise std [deg]
p.sigma_theta_rad = p.sigma_theta_deg * pi/180;   % [rad]

%% Radar blackout bounds
p.t_bo_max    = 0.50;   % maximum blackout duration [s]
p.t_bo0_min   = 1.00;   % earliest blackout onset [s]
p.scnC_t_bo   = 0.50;   % Scenario C blackout duration [s]
p.scnC_t_bo0  = 1.00;   % Scenario C blackout onset [s]

%% Prediction
p.t_p = 0.50;           % prediction horizon [s]

%% Laser
p.R_eng      = 2;       % engagement radius [m]
p.t_recharge = 1.0;     % minimum time between pulses [s]

%% Simulation
p.t_end = 10;           % standardized engagement duration [s]

%% Target operating limits
p.r_min = 800;          % [m]
p.r_max = 2000;         % [m]
p.v_max = 150;          % [m/s]
p.a_max = 20;           % [m/s^2]

%% Scenario B
p.scnB_maneuver_start = 1.0;   % maneuver begins at t = 1 s
end