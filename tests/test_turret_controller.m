function test_turret_controller()
%TEST_TURRET_CONTROLLER  Independent tests for turret_controller.m.
%
%   Run from the folder that holds all project files:   >> test_turret_controller
%
%   Sections
%     1  Design: gains, Routh-Hurwitz, closed-loop poles, gain/phase margins
%     2  Linear loop (real turret_plant, P law): step, ramp lag, parabola residual
%     3  Controller unit tests: init, reset, clamp, sigma_perp, fire/recharge logic
%     4  Static scan: no ground-truth calls in controller/estimator/predictor
%     5  Closed-loop scenarios A, B, C (estimator+predictor+controller+plant)
%     6  Blackout behaviour (Scenario C)
%     7  Angle-limit tests (plant in the loop, targets above/below the limits)
%     8  Figures -> ./figures_controller/*.png
%
%   Two kinds of checks are counted separately:
%     CHECK = structural/design checks that must pass
%     PERF  = performance indicators (a failure is information, not a bug)
%
%   Simulation assumptions (state them in the report):
%     * fixed step dt = 1 ms, t = k*dt computed by multiplication
%     * plant integrated with RK4, controller output held constant over a step (ZOH)
%     * a laser pulse lasts ONE step and is judged at the step where fire_cmd
%       goes 0->1, using the TRUE target position and the turret angle at that step
%     * initial turret angle 0 deg, zero initial speed/current
%     * no voltage/current saturation (none specified)

close all;
fprintf('\n=============== test_turret_controller ===============\n');
chk('__reset__');
t_all = tic;
here   = fileparts(mfilename('fullpath'));
outdir = fullfile(here, 'figures_controller');
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% ------------------------------------------------------------ 0. files
req = {'params.m','turret_plant.m','target_truth.m','target_motion_generator.m', ...
       'target_state_estimator.m','target_predictor.m','turret_controller.m'};
missing = false;
for k = 1:numel(req)
    ok = ~isempty(which(req{k}));
    chk(['file on path: ' req{k}], ok);
    missing = missing || ~ok;
end
if missing
    fprintf('\nMissing files - cd to the project folder and re-run.\n');
    chk('__summary__');
    return;
end
prm = params();
turret_controller();                                  % reset

%% ------------------------------------------------ 1. design / stability
fprintf('\n--- 1. Design: gains, Routh-Hurwitz, poles, margins ---\n');
des = turret_controller('design');
Kv  = prm.G_num / prm.G_den(3);
a3 = prm.G_den(1);  a2 = prm.G_den(2);  a1 = prm.G_den(3);
a0 = prm.G_num * des.Kp;                              % char. eq.: a3 s^3 + a2 s^2 + a1 s + a0
fprintf('  Kv = N*Kt/(Rm*b+N^2*Kt*Ke) = %.6e rad/s/V\n', Kv);
fprintf('  Kp = 1/(Tp*Kv)             = %.4f V/rad   (Tp = %.2f s)\n', des.Kp, prm.t_p);
chk('Kp equals 1/(Tp*Kv)', abs(des.Kp - 1/(prm.t_p*Kv)) < 1e-9*des.Kp);
chk('loop velocity constant Kp*Kv equals 1/Tp', abs(des.Kv_loop - 1/prm.t_p) < 1e-9);

% Routh-Hurwitz for  a3 s^3 + a2 s^2 + a1 s + a0
b1 = (a2*a1 - a3*a0) / a2;
routh_col = [a3; a2; b1; a0];
Kp_max = a2*a1 / (a3*prm.G_num);
fprintf('  Routh first column: [%.4g  %.4g  %.6g  %.6g]\n', routh_col);
fprintf('  Routh limit: Kp < a2*a1/(a3*G_num) = %.2f V/rad  (chosen Kp is %.1f%% of it)\n', ...
        Kp_max, 100*des.Kp/Kp_max);
chk('Routh-Hurwitz: all first-column entries > 0', all(routh_col > 0), ...
    sprintf('s^1 entry = %.6g', b1));
chk('Kp below Routh limit', des.Kp < Kp_max);

cl_poles = roots([a3 a2 a1 a0]);
ol_poles = roots(prm.G_den);
fprintf('  Open-loop poles : '); fprintf('%.3f%+.3fj  ', [real(ol_poles) imag(ol_poles)].'); fprintf('\n');
fprintf('  Closed-loop poles: '); fprintf('%.3f%+.3fj  ', [real(cl_poles) imag(cl_poles)].'); fprintf('\n');
chk('all closed-loop poles in open left half plane', all(real(cl_poles) < 0));
[~, ir] = min(abs(imag(cl_poles)));
p_dom = real(cl_poles(ir));
cpl = cl_poles(abs(imag(cl_poles)) > 1e-6);
if ~isempty(cpl)
    wn_cl = abs(cpl(1));  zeta_cl = -real(cpl(1))/wn_cl;
    fprintf('  Dominant real pole %.3f (tau = %.3f s); resonant pair wn = %.2f rad/s, zeta = %.3f\n', ...
            p_dom, -1/p_dom, wn_cl, zeta_cl);
end

Lfun = @(w) (des.Kp*prm.G_num) ./ ((1i*w) .* polyval(prm.G_den(1:3), 1i*w));
wgc = fzero(@(w) abs(Lfun(w)) - 1, [0.1 20]);
PM  = 180 + angle(Lfun(wgc))*180/pi;
wpc = sqrt(prm.G_den(3)/prm.G_den(1));               % where Im{L}=0 (unique, see report)
GM  = 1/abs(Lfun(wpc));
fprintf('  Gain crossover  wgc = %.3f rad/s, phase margin = %.2f deg\n', wgc, PM);
fprintf('  Phase crossover wpc = %.3f rad/s, gain margin  = %.3f (%.2f dB)\n', wpc, GM, 20*log10(GM));
chk('phase at wpc is -180 deg', abs(abs(angle(Lfun(wpc))) - pi) < 1e-6);
chk('gain margin from margins = Kp_max/Kp (Routh consistency)', abs(GM - Kp_max/des.Kp) < 1e-6*GM);
chk('phase margin > 60 deg', PM > 60, sprintf('PM = %.1f deg', PM));
chk('gain margin > 6 dB',    20*log10(GM) > 6, sprintf('GM = %.1f dB', 20*log10(GM)));

wv = logspace(-1, 3, 4000);  Lw = Lfun(wv);  Tw = Lw ./ (1 + Lw);
ibw = find(abs(Tw) < 1/sqrt(2), 1);  wbw = wv(ibw);
[Tpk, ipk] = max(abs(Tw));
fprintf('  Closed-loop -3 dB bandwidth = %.2f rad/s; resonance peak |T| = %.2f dB at %.1f rad/s\n', ...
        wbw, 20*log10(Tpk), wv(ipk));

if exist('margin','file') && exist('tf','file')
    try
        Gp = tf(prm.G_num, prm.G_den);
        [Gm_t, Pm_t] = margin(des.Kp*Gp);
        fprintf('  Control System Toolbox margin(): GM = %.3f, PM = %.2f deg\n', Gm_t, Pm_t);
        chk('toolbox margin() agrees with manual margins', ...
            abs(Gm_t-GM)/GM < 1e-3 && abs(Pm_t-PM) < 0.1);
    catch err
        fprintf('  (toolbox cross-check skipped: %s)\n', err.message);
    end
else
    fprintf('  (Control System Toolbox not found - manual margins only)\n');
end

%% ------------------------------------ 2. linear loop, real turret_plant
fprintf('\n--- 2. Linear loop with the real turret_plant (P law, no estimator) ---\n');
dtl = 1e-3;
% step
th_step = 0.3;
S1 = sim_linear_loop(@(t) th_step*ones(size(t)), des.Kp, dtl, 10);
m  = step_metrics(S1.t, S1.th, th_step);
fprintf('  Step 0.3 rad: final = %.5f, overshoot = %.2f %%, 10-90%% rise = %.3f s, 2%% settle = %.3f s\n', ...
        m.final, m.overshoot_pct, m.rise, m.settle);
fprintf('                peak V = %.1f V, peak i = %.1f A (no saturation modelled)\n', ...
        max(abs(S1.V)), max(abs(S1.i)));
chk('step: zero steady-state error (type-1)', abs(m.final - th_step) < 1e-4*th_step);
chk('step: overshoot < 1 %', m.overshoot_pct < 1);
chk('step: 2% settling time < 2.5 s', m.settle < 2.5, sprintf('%.3f s', m.settle));
chk('step: all signals finite', all(isfinite(S1.th)) && all(isfinite(S1.V)) && all(isfinite(S1.i)));

% ramp: steady-state lag = 1/Kv_loop = Tp
rate = 0.04;
S2 = sim_linear_loop(@(t) rate*t, des.Kp, dtl, 10);
late = S2.t >= 8;
lag  = mean((S2.ref(late) - S2.th(late))) / rate;
fprintf('  Ramp %.3f rad/s: steady-state lag = %.4f s   (theory 1/(Kp*Kv) = %.4f s)\n', rate, lag, 1/des.Kv_loop);
chk('ramp: steady-state lag equals 1/(Kp*Kv) = Tp', abs(lag - 1/des.Kv_loop) < 0.01, sprintf('%.4f s', lag));

% parabola: pointing error w.r.t. the CURRENT target, reference = target(t+Tp)
alpha = 0.02;  w0 = 0.03;  th0q = 0.3;
tgt = @(t) th0q + w0*t + 0.5*alpha*t.^2;
S3  = sim_linear_loop(@(t) tgt(t + prm.t_p), des.Kp, dtl, 10);
c1 = a1/a0;  c2 = a2/a0;
c_res = prm.t_p^2/2 - c1*prm.t_p + (c1^2 - c2);
e_pt  = mean(S3.th(late) - tgt(S3.t(late)));
fprintf('  Constant angular acceleration %.3f rad/s^2: pointing error = %.3e rad (theory c*alpha = %.3e, c = %.4f)\n', ...
        alpha, e_pt, c_res*alpha, c_res);
chk('parabola: residual pointing error = c*theta_ddot (no perfect cancellation)', ...
    abs(e_pt - c_res*alpha) < 0.01*c_res*alpha);

% how big could that residual be in the scenarios (quasi-static estimate)?
tt = (0:0.01:prm.t_end)';
for sc = {'A','B'}
    [xs, ys] = target_truth(tt, sc{1});
    ths = atan2(ys, xs);  rs = hypot(xs, ys);
    thd = gradient(ths, 0.01);  thdd = gradient(thd, 0.01);
    fprintf('  Scenario %s (truth, test-only): max|theta_dot| = %.4f rad/s, max|theta_ddot| = %.4f rad/s^2, quasi-static residual max = %.2f m\n', ...
            sc{1}, max(abs(thd)), max(abs(thdd)), max(abs(c_res*thdd.*rs)));
end

%% -------------------------------------------- 3. controller unit tests
fprintf('\n--- 3. Controller unit tests (synthetic inputs) ---\n');

% U1: uninitialised
turret_controller();
[V, f] = turret_controller(0, NaN, NaN, false, 0, 0);
inf0 = turret_controller('info');
chk('no measurement yet: V_m = 0, no fire, estimator not valid', V == 0 && ~f && ~inf0.valid);

% U2: initialisation on a stationary target, turret pointing at it
turret_controller();
r_s = 1000;  th_s = 0.3;
[V, f] = turret_controller(0, r_s, th_s, true, th_s, 0);
i2 = turret_controller('info');
chk('init: finite V_m, essentially zero (reference = current angle)', isfinite(V) && abs(V) < 1e-6, sprintf('V = %.3g', V));
chk('init: fire_cmd high when turret points at target (gate open)', f == true, ...
    sprintf('d_hat = %.3g m, sigma_perp = %.3f m', i2.d_hat, i2.sigma_perp));
chk('sigma_perp = n''*Pxy*n = r*sigma_theta when n is the cross-range direction', ...
    abs(i2.sigma_perp - r_s*prm.sigma_theta_rad) < 1e-9*r_s*prm.sigma_theta_rad, ...
    sprintf('%.6f vs %.6f m', i2.sigma_perp, r_s*prm.sigma_theta_rad));

% U3: sigma_perp when the beam normal points along the range direction
turret_controller();
[V, f] = turret_controller(0, r_s, th_s, true, th_s + pi/2, 0);
i3 = turret_controller('info');
chk('sigma_perp = sigma_r when beam normal is along range', abs(i3.sigma_perp - prm.sigma_r) < 1e-9, ...
    sprintf('%.6f m', i3.sigma_perp));
chk('beam 90 deg off target: no fire, finite V', ~f && isfinite(V));

% U4: reference clamp (turret limits)
turret_controller();
[V, ~] = turret_controller(0, 1000, 89*pi/180, true, 0.2, 0);
i4 = turret_controller('info');
chk('reference clamped at +85 deg', abs(i4.theta_ref - prm.theta_max_rad) < 1e-12, ...
    sprintf('raw %.2f deg -> %.2f deg', i4.theta_ref_raw*180/pi, i4.theta_ref*180/pi));
chk('V_m uses the clamped reference', abs(V - des.Kp*(prm.theta_max_rad - 0.2)) < 1e-9*abs(V));
turret_controller();
[~, ~] = turret_controller(0, 1000, -25*pi/180, true, 0.0, 0);
i5 = turret_controller('info');
chk('reference clamped at -10 deg', abs(i5.theta_ref - prm.theta_min_rad) < 1e-12, ...
    sprintf('raw %.2f deg -> %.2f deg', i5.theta_ref_raw*180/pi, i5.theta_ref*180/pi));

% U5: reset when time moves backwards (also clears recharge timer)
turret_controller();
for k = 0:20, turret_controller(k*0.05, r_s, th_s, true, th_s, 0); end
i6a = turret_controller('info');
[~, f] = turret_controller(0, r_s, th_s, true, th_s, 0);        % time goes backwards
i6b = turret_controller('info');
chk('21 radar updates counted before the reset', i6a.est_n_updates == 21);
chk('time moving backwards restarts estimator (n_updates = 1, t_init = 0)', ...
    i6b.est_n_updates == 1 && i6b.est_t_init == 0);
chk('reset also clears the recharge timer (pulse allowed at new t = 0)', f == true);

% U6: laser pulse transitions and recharge timing (stationary target, turret on it)
[~, Fs] = run_stationary(45, 0.05, r_s, th_s, th_s);
tpl = (0:44)*0.05;  tpl = tpl(Fs);
fprintf('  On-target stationary case, pulse times: %s s\n', mat2str(tpl, 4));
chk('pulse train: first pulse at t = 0', ~isempty(tpl) && tpl(1) == 0);
chk('pulse spacing >= 1.0 s (recharge)', numel(tpl) < 2 || all(diff(tpl) >= prm.t_recharge));
chk('pulse train: 3 pulses at 0, 1.05, 2.10 s (1.0 s itself blocked by round-off guard)', ...
    numel(tpl) == 3 && all(abs(tpl - [0 1.05 2.10]) < 1e-9));
chk('fire_cmd never high on two consecutive steps', ~any(Fs(1:end-1) & Fs(2:end)));
[~, Fs2] = run_stationary(21, 0.05, r_s, th_s, th_s + 5*pi/180);
chk('beam 5 deg off target: gate stays closed, no pulse', ~any(Fs2));

% U7: gate closed by uncertainty even though the estimate is on the beam
turret_controller();
turret_controller(0, r_s, th_s, true, th_s + 5*pi/180, 0);       % init, off boresight -> no pulse
[~, f] = turret_controller(1.2, NaN, NaN, false, th_s, 0);       % 1.2 s blind, turret on estimate
i7 = turret_controller('info');
chk('1.2 s without measurements: estimate on beam but sigma_perp > R_eng -> no fire', ...
    ~f && i7.d_hat < prm.R_eng && i7.sigma_perp > prm.R_eng, ...
    sprintf('d_hat = %.2f m, sigma_perp = %.1f m', i7.d_hat, i7.sigma_perp));

% U8: determinism
[V1, F1] = run_stationary(30, 0.05, r_s, th_s, 0.25);
[V2, F2] = run_stationary(30, 0.05, r_s, th_s, 0.25);
chk('identical inputs after reset give identical outputs', isequal(V1, V2) && isequal(F1, F2));

%% ------------------------------------- 4. no ground truth in the chain
fprintf('\n--- 4. Static scan for ground-truth use ---\n');
banned = '(target_truth|scenario_definition|target_motion_generator)';
files  = {'turret_controller.m','target_state_estimator.m','target_predictor.m'};
for k = 1:numel(files)
    txt   = fileread(which(files{k}));
    lines = regexp(txt, '\r?\n', 'split');
    hits  = 0;
    for j = 1:numel(lines)
        code = regexprep(lines{j}, '%.*$', '');      % strip comments
        if ~isempty(regexp(code, banned, 'once')), hits = hits + 1; end
    end
    chk(['no ground-truth call in code of ' files{k}], hits == 0, sprintf('%d offending lines', hits));
end

%% ------------------------------------------------ 5. closed-loop A, B, C
fprintf('\n--- 5. Closed-loop scenarios (dt = 1 ms, RK4 plant, ZOH controller) ---\n');
dt = 1e-3;  th_init = 0;
scn = {'A','B','C'};
Ls = cell(1,3);  Ss = cell(1,3);
for s = 1:3
    tic;
    Ls{s} = run_scenario(scn{s}, dt, th_init);
    Ss{s} = summarize(Ls{s}, prm);
    fprintf('\n  Scenario %s  (%.1f s wall time)\n', scn{s}, toc);
    print_summary(Ss{s}, prm);
    L = Ls{s};  S = Ss{s};
    tag = ['Scn ' scn{s} ': '];
    chk([tag 'V_m, turret angle, current all finite'], all(isfinite(L.V)) && all(isfinite(L.th)) && all(isfinite(L.i)));
    chk([tag 'turret stays inside [-10, 85] deg'], ...
        min(L.th) >= prm.theta_min_rad - 1e-6 && max(L.th) <= prm.theta_max_rad + 1e-6, ...
        sprintf('%.2f .. %.2f deg', min(L.th)*180/pi, max(L.th)*180/pi));
    chk([tag 'first command is Kp*(ref - theta0) (pure P law, no kick)'], ...
        abs(L.V(1) - des.Kp*(L.th_ref(1) - th_init)) < 1e-9*abs(L.V(1)));
    chk([tag 'estimator never reset (t_init = 0 throughout)'], all(L.tinit == 0));
    chk([tag 'no pulse transition ignored by the recharge rule'], L.ignored == 0);
    chk([tag 'pulse spacing >= t_recharge'], numel(L.pulse_t) < 2 || all(diff(L.pulse_t) >= prm.t_recharge));
    chk([tag 'fire_cmd never high on consecutive steps'], ~any(L.fire(1:end-1) & L.fire(2:end)));
    if s < 3
        chk([tag 'all 201 radar samples used by estimator'], L.nupd(end) == 201, sprintf('%d', L.nupd(end)));
    else
        chk([tag '191 radar samples used (10 lost in blackout)'], L.nupd(end) == 191, sprintf('%d', L.nupd(end)));
    end
    chk([tag 'at least one pulse issued'], S.npulse >= 1, 'PERF-style', 'perf');
    chk([tag 'at least one successful engagement (<= R_eng at firing instant)'], S.nhit >= 1, ...
        sprintf('%d hit(s) of %d pulse(s)', S.nhit, S.npulse), 'perf');
end

%% ----------------------------------------------------- 6. blackout (C)
fprintf('\n--- 6. Blackout behaviour (Scenario C, samples 1.00-1.45 s lost) ---\n');
LC = Ls{3};  LB = Ls{2};
ix = @(tq) round(tq/dt) + 1;
n_before = LC.nupd(ix(0.95));  n_end = LC.nupd(ix(1.499));  n_back = LC.nupd(ix(1.5));
chk('estimator count frozen during blackout', n_before == n_end, sprintf('%d -> %d', n_before, n_end));
chk('estimator resumes updating at 1.50 s', n_back == n_end + 1, sprintf('%d -> %d', n_end, n_back));
chk('sigma_perp grows during blackout', LC.sig(ix(1.499)) > LC.sig(ix(0.95)), ...
    sprintf('%.2f m -> %.2f m', LC.sig(ix(0.95)), LC.sig(ix(1.499))));
bo = LC.t >= 1.0 & LC.t < 1.5;
chk('blackout: V_m finite and continuous (no reset/glitch)', all(isfinite(LC.V(bo))) && max(abs(diff(LC.V(bo)))) < 50, ...
    sprintf('max step in V = %.2f V', max(abs(diff(LC.V(bo))))));
fprintf('  Gate open during blackout: %.1f %% of steps (C) vs %.1f %% (B, same motion, no blackout)\n', ...
        100*mean(LC.gate(bo)), 100*mean(LB.gate(bo)));
jump = abs(LC.th_ref(ix(1.5)) - LC.th_ref(ix(1.499)));
jumpB = abs(LB.th_ref(ix(1.5)) - LB.th_ref(ix(1.499)));
fprintf('  Reference jump at measurement recovery: C = %.3e rad, B = %.3e rad\n', jump, jumpB);
w1 = LC.t >= 1.0 & LC.t <= 2.0 & isfinite(LC.pred_err);
fprintf('  Predicted-elevation RMS error over 1-2 s: C = %.3e rad, B = %.3e rad\n', ...
        sqrt(mean(LC.pred_err(w1).^2)), sqrt(mean(LB.pred_err(w1).^2)));
fprintf('  Pointing miss (true) max over 1-2 s: C = %.2f m, B = %.2f m\n', ...
        max(LC.miss(w1)), max(LB.miss(w1)));

%% ------------------------------------------------- 7. angle-limit tests
fprintf('\n--- 7. Angle-limit tests (plant in the loop) ---\n');
Lu = run_synthetic(1000, 89*pi/180, 0, 6, dt);
chk('target above limit: turret never exceeds +85 deg', max(Lu.th) <= prm.theta_max_rad + 1e-4, ...
    sprintf('max = %.4f deg', max(Lu.th)*180/pi));
chk('target above limit: turret reaches the 85 deg stop (within 0.05 deg)', ...
    abs(Lu.th(end) - prm.theta_max_rad) < 0.05*pi/180, sprintf('final = %.4f deg', Lu.th(end)*180/pi));
Ll = run_synthetic(1000, -25*pi/180, 0, 6, dt);
chk('target below limit: turret never goes below -10 deg', min(Ll.th) >= prm.theta_min_rad - 1e-4, ...
    sprintf('min = %.4f deg', min(Ll.th)*180/pi));
chk('target below limit: turret reaches -10 deg (within 0.05 deg)', ...
    abs(Ll.th(end) - prm.theta_min_rad) < 0.05*pi/180, sprintf('final = %.4f deg', Ll.th(end)*180/pi));
chk('limit tests: finite signals', all(isfinite(Lu.V)) && all(isfinite(Ll.V)));
fprintf('  Initial slew to 85 deg: peak V = %.0f V, peak i = %.0f A (no saturation modelled)\n', ...
        max(abs(Lu.V)), max(abs(Lu.i)));

%% ---------------------------------------------------------- 8. figures
fprintf('\n--- 8. Figures ---\n');
col = {'b','r','m'};
f1 = figure('Name','Elevation tracking','Position',[50 50 900 800]);
for s = 1:3
    L = Ls{s};
    subplot(3,1,s);
    plot(L.t, L.th_true*180/pi, 'k', 'LineWidth', 1.2); hold on;
    plot(L.t, L.th_ref*180/pi, 'b--');
    plot(L.t, L.th*180/pi, 'r');
    if ~isempty(L.pulse_t)
        h = L.pulse_hit;
        plot(L.pulse_t(h),  L.pulse_th(h)*180/pi,  'g^', 'MarkerFaceColor','g', 'MarkerSize', 8);
        plot(L.pulse_t(~h), L.pulse_th(~h)*180/pi, 'rx', 'MarkerSize', 9, 'LineWidth', 2);
    end
    if s == 3, patch_blackout(); end
    grid on; ylabel('elevation [deg]'); title(sprintf('Scenario %s', scn{s}));
    if s == 1, legend('true now','reference (pred. t+T_p)','turret','hit','miss','Location','best'); end
end
xlabel('t [s]');  save_fig(f1, outdir, 'fig1_elevation_tracking');

f2 = figure('Name','Pointing miss distance','Position',[60 60 900 800]);
for s = 1:3
    L = Ls{s};
    subplot(3,1,s);
    plot(L.t, L.miss, col{s}); hold on;
    plot([0 prm.t_end], [prm.R_eng prm.R_eng], 'k--');
    if ~isempty(L.pulse_t), plot(L.pulse_t, min(L.pulse_miss, 40), 'ko', 'MarkerFaceColor','y'); end
    if s == 3, patch_blackout(); end
    set(gca, 'YLim', [0 40]); grid on; ylabel('miss [m]');
    title(sprintf('Scenario %s: beam-to-target distance (true), pulses = yellow dots, R_{eng} dashed', scn{s}));
end
xlabel('t [s]');  save_fig(f2, outdir, 'fig2_pointing_miss');

f3 = figure('Name','Control effort','Position',[70 70 900 600]);
subplot(2,1,1); hold on; grid on;
for s = 1:3, plot(Ls{s}.t, Ls{s}.V, col{s}); end
ylabel('V_m [V]'); legend(scn, 'Location','best'); title('Motor voltage (no saturation modelled)');
subplot(2,1,2); hold on; grid on;
for s = 1:3, plot(Ls{s}.t, Ls{s}.i, col{s}); end
ylabel('i [A]'); xlabel('t [s]');
save_fig(f3, outdir, 'fig3_control_effort');

f4 = figure('Name','Blackout zoom (C)','Position',[80 80 900 700]);
zz = LC.t >= 0.8 & LC.t <= 2.2;
subplot(3,1,1); hold on; grid on;
plot(LC.t(zz), LC.th_ref_raw(zz)*180/pi, 'b');
plot(LC.t(zz), (LC.th_true(zz))*180/pi, 'k');
plot(LC.t(zz), LB.th_ref_raw(zz)*180/pi, 'g--');
patch_blackout(); ylabel('[deg]'); legend('C: predicted ref','true now','B: predicted ref','Location','best');
title('Predicted elevation reference across the blackout');
subplot(3,1,2); hold on; grid on;
semilogy(LC.t(zz), LC.sig(zz), 'r'); set(gca,'YScale','log');
plot([0.8 2.2], [prm.R_eng prm.R_eng], 'k--'); patch_blackout();
ylabel('\sigma_{perp} [m]'); title('Cross-beam position uncertainty (n^T P_{xy} n)^{1/2}');
subplot(3,1,3); hold on; grid on;
plot(LC.t(zz), LC.d_hat(zz) + LC.sig(zz), 'm');
plot([0.8 2.2], [prm.R_eng prm.R_eng], 'k--'); patch_blackout();
set(gca,'YLim',[0 12]); ylabel('d_{hat}+\sigma_{perp} [m]'); xlabel('t [s]');
title('Fire gate: open when below the dashed R_{eng} line');
save_fig(f4, outdir, 'fig4_blackout_zoom');

f5 = figure('Name','Linear loop analysis','Position',[90 90 1000 800]);
subplot(2,2,1);
semilogx(wv, 20*log10(abs(Lw)), 'b'); hold on; grid on; plot([wv(1) wv(end)],[0 0],'k:');
plot(wgc, 0, 'ro'); plot(wpc, -20*log10(GM), 'rs');
ylabel('|L| [dB]'); title('Loop gain L = K_p G_p'); xlim([wv(1) wv(end)]);
subplot(2,2,3);
semilogx(wv, unwrap(angle(Lw))*180/pi, 'b'); hold on; grid on; plot([wv(1) wv(end)],[-180 -180],'k:');
plot(wgc, -180+PM, 'ro'); xlabel('\omega [rad/s]'); ylabel('phase [deg]'); xlim([wv(1) wv(end)]);
subplot(2,2,2);
plot(S1.t, S1.th, 'b'); hold on; plot(S1.t, th_step*ones(size(S1.t)), 'k--'); grid on;
xlim([0 4]); xlabel('t [s]'); ylabel('\theta [rad]'); title('Step response (real turret\_plant)');
subplot(2,2,4);
plot(S2.t, S2.ref - S2.th, 'r'); hold on; plot([0 10],[rate/des.Kv_loop rate/des.Kv_loop],'k--'); grid on;
xlabel('t [s]'); ylabel('ramp error [rad]'); title('Ramp error -> rate/(K_pK_v) = rate\cdotT_p');
save_fig(f5, outdir, 'fig5_linear_loop');
fprintf('  PNGs written to %s\n', outdir);

%% ------------------------------------------------------------- summary
fprintf('\nTotal wall time: %.1f s\n', toc(t_all));
chk('__summary__');
end

% =====================================================================
%                              LOCAL FUNCTIONS
% =====================================================================
function chk(name, cond, detail, cat)
% Records and prints one check.  chk('__reset__') / chk('__summary__') are service calls.
persistent store
if isempty(store), store = struct('name', {}, 'ok', {}, 'cat', {}); end
if strcmp(name, '__reset__'), store = struct('name', {}, 'ok', {}, 'cat', {}); return; end
if strcmp(name, '__summary__')
    if isempty(store), fprintf('No checks recorded.\n'); return; end
    isperf = strcmp({store.cat}, 'perf');  ok = [store.ok];
    fprintf('\n================ SUMMARY ================\n');
    fprintf('CHECK (structural/design): %d / %d passed\n', sum(ok & ~isperf), sum(~isperf));
    fprintf('PERF  (performance flags): %d / %d passed\n', sum(ok & isperf),  sum(isperf));
    bad = find(~ok);
    if isempty(bad)
        fprintf('No failures.\n');
    else
        fprintf('Failures:\n');
        for k = bad, fprintf('  - [%s] %s\n', upper(store(k).cat), store(k).name); end
    end
    return;
end
if nargin < 3, detail = ''; end
if nargin < 4, cat = 'check'; end
ok = ~isempty(cond) && all(cond(:));
k = numel(store) + 1;
store(k).name = name;  store(k).ok = ok;  store(k).cat = cat;
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
if strcmp(cat, 'perf'), tag = ['PERF-' tag]; end
if isempty(detail)
    fprintf('  [%s] %s\n', tag, name);
else
    fprintf('  [%s] %s  (%s)\n', tag, name, detail);
end
end

function xd = plant_deriv(t, x, V)
[w, a, di] = turret_plant(t, V, x(1), x(2), x(3));
xd = [w; a; di];
end

function x = rk4_step(x, V, t, dt)
k1 = plant_deriv(t,        x,             V);
k2 = plant_deriv(t + dt/2, x + dt/2*k1,   V);
k3 = plant_deriv(t + dt/2, x + dt/2*k2,   V);
k4 = plant_deriv(t + dt,   x + dt*k3,     V);
x  = x + dt/6*(k1 + 2*k2 + 2*k3 + k4);
end

function S = sim_linear_loop(refFun, Kp, dt, Tend)
% P-law loop around the real turret_plant, starting at rest at theta = 0.
N = round(Tend/dt);  x = [0; 0; 0];
S.t = zeros(N+1,1); S.th = S.t; S.V = S.t; S.i = S.t; S.ref = S.t;
for k = 0:N
    t = k*dt;  r = refFun(t);  V = Kp*(r - x(1));
    S.t(k+1) = t; S.th(k+1) = x(1); S.i(k+1) = x(3); S.V(k+1) = V; S.ref(k+1) = r;
    if k < N, x = rk4_step(x, V, t, dt); end
end
end

function m = step_metrics(t, y, yf)
m.final = y(end);
m.overshoot_pct = max(0, (max(y) - yf)/yf*100);
i10 = find(y >= 0.1*yf, 1);  i90 = find(y >= 0.9*yf, 1);
m.rise = t(i90) - t(i10);
idx = find(abs(y - yf) > 0.02*abs(yf), 1, 'last');
if isempty(idx), m.settle = 0; else, m.settle = t(min(idx+1, numel(t))); end
end

function d = ray_miss(x, y, theta)
% Distance from point (x,y) to the beam ray (origin, direction theta).
along = x.*cos(theta) + y.*sin(theta);
d = abs(-x.*sin(theta) + y.*cos(theta));
back = along <= 0;
d(back) = hypot(x(back), y(back));
end

function [Vs, Fs] = run_stationary(nsteps, dt, r_s, th_s, th_turret)
% Feeds the controller noise-free measurements of a stationary target.
turret_controller();
Vs = zeros(nsteps,1);  Fs = false(nsteps,1);
for k = 1:nsteps
    [Vs(k), f] = turret_controller((k-1)*dt, r_s, th_s, true, th_turret, 0);
    Fs(k) = f;
end
end

function L = run_synthetic(r_s, th_s, th0, Tend, dt)
% Plant in the loop, stationary synthetic target measured at 20 Hz.
turret_controller();
N = round(Tend/dt);  per = round(0.05/dt);  x = [th0; 0; 0];
L.t = zeros(N+1,1); L.th = L.t; L.V = L.t; L.i = L.t;
for k = 0:N
    t = k*dt;  valid = (mod(k, per) == 0);
    if valid, rm = r_s; tm = th_s; else, rm = NaN; tm = NaN; end
    V = turret_controller(t, rm, tm, valid, x(1), x(2));
    L.t(k+1) = t; L.th(k+1) = x(1); L.V(k+1) = V; L.i(k+1) = x(3);
    if k < N, x = rk4_step(x, V, t, dt); end
end
end

function L = run_scenario(scn, dt, th0)
prm = params();
turret_controller();
N = round(prm.t_end/dt);  n = N + 1;
L.t = zeros(n,1);  L.th = zeros(n,1);  L.w = zeros(n,1);  L.i = zeros(n,1);  L.V = zeros(n,1);
L.th_true = zeros(n,1);  L.r_true = zeros(n,1);  L.th_ref = zeros(n,1);  L.th_ref_raw = nan(n,1);
L.miss = zeros(n,1);  L.d_hat = nan(n,1);  L.sig = nan(n,1);
L.gate = false(n,1);  L.fire = false(n,1);  L.pred_err = nan(n,1);
L.nupd = zeros(n,1);  L.tinit = nan(n,1);
L.pulse_t = [];  L.pulse_hit = false(0,1);  L.pulse_miss = [];  L.pulse_th = [];  L.ignored = 0;
x = [th0; 0; 0];  prev_fire = false;  last_pulse = -Inf;
for k = 0:N
    t = k*dt;
    [r_m, th_m, v_m] = target_motion_generator(t, scn);
    [V, fire] = turret_controller(t, r_m, th_m, v_m, x(1), x(2));
    inf_ = turret_controller('info');
    [xt, yt] = target_truth(t, scn);                      % TEST-ONLY evaluation
    miss = ray_miss(xt, yt, x(1));
    fire = logical(fire);
    if fire && ~prev_fire                                  % 0 -> 1 transition
        if (t - last_pulse) >= prm.t_recharge
            last_pulse = t;
            L.pulse_t(end+1,1) = t;
            L.pulse_hit(end+1,1) = (miss <= prm.R_eng);
            L.pulse_miss(end+1,1) = miss;
            L.pulse_th(end+1,1) = x(1);
        else
            L.ignored = L.ignored + 1;
        end
    end
    prev_fire = fire;
    j = k + 1;
    L.t(j) = t;  L.th(j) = x(1);  L.w(j) = x(2);  L.i(j) = x(3);  L.V(j) = V;
    L.th_true(j) = atan2(yt, xt);  L.r_true(j) = hypot(xt, yt);  L.miss(j) = miss;
    L.fire(j) = fire;  L.nupd(j) = inf_.est_n_updates;  L.tinit(j) = inf_.est_t_init;
    if inf_.valid
        [xf, yf] = target_truth(t + prm.t_p, scn);
        L.th_ref(j) = inf_.theta_ref;  L.th_ref_raw(j) = inf_.theta_ref_raw;
        L.d_hat(j) = inf_.d_hat;  L.sig(j) = inf_.sigma_perp;  L.gate(j) = inf_.gate_open;
        L.pred_err(j) = inf_.theta_ref_raw - atan2(yf, xf);
    end
    if k < N, x = rk4_step(x, V, t, dt); end
end
L.dt = dt;  L.scn = scn;
end

function S = summarize(L, prm)
t = L.t;  w3 = t >= 3;
S.rms_miss  = sqrt(mean(L.miss(w3).^2));
S.max_miss  = max(L.miss(w3));
S.frac_in   = mean(L.miss(w3) <= prm.R_eng);
ia = find(L.miss <= prm.R_eng, 1);
if isempty(ia), S.t_acq = NaN; else, S.t_acq = t(ia); end
S.npulse    = numel(L.pulse_t);
S.nhit      = sum(L.pulse_hit);
S.first_hit = NaN;
if any(L.pulse_hit), S.first_hit = L.pulse_t(find(L.pulse_hit, 1)); end
S.Vpk = max(abs(L.V));  S.Vpk_1s = max(abs(L.V(t >= 1)));  S.Vrms = sqrt(mean(L.V.^2));
S.Ipk = max(abs(L.i));  S.Ipk_1s = max(abs(L.i(t >= 1)));
S.th_min = min(L.th)*180/pi;  S.th_max = max(L.th)*180/pi;
ok = isfinite(L.pred_err) & t >= 1;
S.pred_rms_rad = sqrt(mean(L.pred_err(ok).^2));
S.pred_rms_m   = sqrt(mean((L.pred_err(ok).*L.r_true(ok)).^2));
S.gate_frac    = mean(L.gate);
end

function print_summary(S, prm)
fprintf('    first time beam within R_eng of target ......... %.3f s\n', S.t_acq);
fprintf('    pointing miss, t >= 3 s: RMS %.2f m, max %.2f m, within R_eng %.1f %% of time\n', ...
        S.rms_miss, S.max_miss, 100*S.frac_in);
fprintf('    pulses: %d issued, %d hit; first hit at %.3f s (NaN = none)\n', S.npulse, S.nhit, S.first_hit);
fprintf('    predicted-elevation error (t >= 1 s): RMS %.3e rad = %.2f m\n', S.pred_rms_rad, S.pred_rms_m);
fprintf('    control effort: |V| peak %.0f V (%.1f V after 1 s), RMS %.1f V; |i| peak %.0f A (%.1f A after 1 s)\n', ...
        S.Vpk, S.Vpk_1s, S.Vrms, S.Ipk, S.Ipk_1s);
fprintf('    turret angle range: %.2f .. %.2f deg (limits %.0f .. %.0f)\n', S.th_min, S.th_max, ...
        prm.theta_min_deg, prm.theta_max_deg);
fprintf('    fire gate open %.1f %% of steps\n', 100*S.gate_frac);
end

function patch_blackout()
yl = get(gca, 'YLim');
p = patch([1 1.5 1.5 1], [yl(1) yl(1) yl(2) yl(2)], [0.8 0.8 0.8], 'FaceAlpha', 0.35, 'EdgeColor', 'none');
try
    set(get(get(p,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
catch
end
set(gca, 'YLim', yl);
end

function save_fig(fig, outdir, name)
try
    saveas(fig, fullfile(outdir, [name '.png']));
catch err
    fprintf('  (could not save %s: %s)\n', name, err.message);
end
end
