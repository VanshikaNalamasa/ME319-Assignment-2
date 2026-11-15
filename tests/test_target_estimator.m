% test_target_estimator.m
% Tests target_state_estimator.m and target_predictor.m.
% Measurements come ONLY from target_motion_generator. Ground truth
% (target_truth) is used ONLY here, for scoring, never passed to the estimator.
% Requires: params.m, target_truth.m, target_motion_generator.m,
%           target_state_estimator.m, target_predictor.m in the same folder.
clear; clc; close all;

P      = params();
scn    = {'A','B','C'};
dt_sim = 0.01;                                   % controller-style call step [s]
t_sim  = (0:round(P.t_end/dt_sim))' * dt_sim;
nS     = numel(t_sim);
checks = cell(0,3);
fprintf('=== test_target_estimator ===\n');
fprintf('Estimator called every %.0f ms; radar at 20 Hz.\n\n', dt_sim*1000);

%% 1) Measurement sequences (generator is the only data source)
for s = 1:3
    M(s).r = nan(nS,1);  M(s).th = nan(nS,1);  M(s).valid = false(nS,1);
    for j = 1:nS
        [M(s).r(j), M(s).th(j), M(s).valid(j)] = target_motion_generator(t_sim(j), scn{s});
    end
end

%% 2) Run estimator + predictor on each scenario
est_end = cell(1,3);
for s = 1:3
    [h_tmp, e_tmp] = run_filter(t_sim, M(s).r, M(s).th, M(s).valid, NaN, true);
    H(s) = h_tmp;  est_end{s} = e_tmp;
end

%% 3) Ground truth for scoring (test script only)
for s = 1:3
    [xt,yt,vxt,vyt,axt,ayt] = target_truth(t_sim, scn{s});
    [xf,yf]                 = target_truth(t_sim + P.t_p, scn{s});
    T(s).x = xt;  T(s).y = yt;  T(s).vx = vxt;  T(s).vy = vyt;  T(s).ax = axt;  T(s).ay = ayt;
    T(s).xf = xf; T(s).yf = yf;
    T(s).pos_err  = hypot(H(s).s(:,1)-xt,  H(s).s(:,2)-yt);
    T(s).vel_err  = hypot(H(s).s(:,3)-vxt, H(s).s(:,4)-vyt);
    T(s).acc_err  = hypot(H(s).s(:,5)-axt, H(s).s(:,6)-ayt);
    T(s).pred_err = hypot(H(s).px-xf, H(s).py-yf);
    T(s).dth      = wrap_pi(H(s).pth - atan2(yf,xf));
    T(s).miss     = hypot(xf,yf) .* abs(sin(T(s).dth));   % miss distance at true future range
end

%% 4) HARD PASS/FAIL CHECKS
fprintf('---- HARD CHECKS (pass/fail) ----\n');
for s = 1:3
    h = H(s);  nm = scn{s};
    r0 = M(s).r(1);  th0 = M(s).th(1);
    s0_ref = [r0*cos(th0); r0*sin(th0); 0; 0; 0; 0];

    checks = add_check(checks, ['Initialises at first valid measurement, ' nm], ...
        M(s).valid(1) && h.init(1) && h.nupd(1)==1 && abs(h.tinit(1)-t_sim(1))<1e-12, ...
        sprintf('t_init = %.2f s', h.tinit(1)));
    checks = add_check(checks, ['Initial state = polar->Cartesian of 1st measurement, ' nm], ...
        max(abs(h.s(1,:)' - s0_ref)) < 1e-9, '');
    checks = add_check(checks, ['Estimates and covariance finite at every step, ' nm], ...
        all(isfinite(h.s(:))) && all(h.Pfin), '');
    checks = add_check(checks, ['Covariance symmetric and positive semi-definite, ' nm], ...
        all(h.P_ok), '');
    checks = add_check(checks, ['Estimator time equals call time at every step, ' nm], ...
        max(abs(h.t_est - t_sim)) < 1e-12, '');
    checks = add_check(checks, ['Each valid measurement used exactly once, ' nm], ...
        h.nupd(end) == nnz(M(s).valid), ...
        sprintf('%d updates, %d valid measurements', h.nupd(end), nnz(M(s).valid)));
    checks = add_check(checks, ['Update flag true only on valid-measurement steps, ' nm], ...
        isequal(h.updated, M(s).valid), '');
    checks = add_check(checks, ['Predictor valid and finite elevation at every step, ' nm], ...
        all(h.pvalid) && all(isfinite(h.pth)), '');
    checks = add_check(checks, ['Predictor target time = t + Tp, ' nm], ...
        max(abs(h.ptt - (t_sim + P.t_p))) < 1e-12, '');
    xp_ref = h.s(:,1) + P.t_p*h.s(:,3) + 0.5*P.t_p^2*h.s(:,5);
    yp_ref = h.s(:,2) + P.t_p*h.s(:,4) + 0.5*P.t_p^2*h.s(:,6);
    checks = add_check(checks, ['Prediction = CA extrapolation of the estimate, ' nm], ...
        max(abs([h.px-xp_ref; h.py-yp_ref; h.pth-atan2(yp_ref,xp_ref)])) < 1e-9, '');
end

% ---- Scenario C blackout behaviour (blackout times used ONLY by this test) ----
h  = H(3);
t0 = P.scnC_t_bo0;  t1 = t0 + P.scnC_t_bo;
i_before = find(abs(t_sim - (t0 - 0.05)) < 1e-9);          % last sample before loss
i_in     = find(t_sim >= t0 - 1e-9 & t_sim < t1 - 1e-9);   % blackout steps
i_last   = i_in(end);
i_resume = find(abs(t_sim - t1) < 1e-9);                   % first sample after loss
seg = h.s(i_before:i_resume,:);

checks = add_check(checks, 'C: no measurement update during blackout', ...
    all(h.nupd(i_in) == h.nupd(i_before)) && ~any(h.updated(i_in)), ...
    sprintf('n_updates stays %d', h.nupd(i_before)));
checks = add_check(checks, 'C: estimator keeps propagating (time advances, estimate moves)', ...
    all(diff(h.t_est(i_before:i_resume)) > 0) && norm(h.s(i_last,1:2)-h.s(i_before,1:2)) > 0, ...
    sprintf('est. moved %.1f m during blackout', norm(h.s(i_last,1:2)-h.s(i_before,1:2))));
checks = add_check(checks, 'C: position uncertainty grows during blackout', ...
    h.Ppos(i_last) > h.Ppos(i_before), ...
    sprintf('pos sigma %.2f m -> %.2f m', sqrt(h.Ppos(i_before)), sqrt(h.Ppos(i_last))));
checks = add_check(checks, 'C: update resumes with first post-blackout measurement', ...
    h.updated(i_resume) && h.nupd(i_resume) == h.nupd(i_before)+1, '');
checks = add_check(checks, 'C: no reset, estimate finite through blackout', ...
    all(h.tinit(i_before:i_resume) == t_sim(1)) && all(isfinite(seg(:))), '');

% ---- Unit tests of estimator behaviour ----
e = target_state_estimator();
e = target_state_estimator(e, 0.00, NaN, NaN, false);
p_u = target_predictor(e);
checks = add_check(checks, 'Uninitialised: no crash, predictor reports invalid', ...
    ~e.initialized && ~p_u.valid && isnan(p_u.theta_p), '');
e = target_state_estimator(e, 0.05, 1500, 0.4, true);
checks = add_check(checks, 'Initialises when first valid measurement arrives', ...
    e.initialized && e.n_updates==1 && abs(e.s(1)-1500*cos(0.4))<1e-9 && abs(e.s(2)-1500*sin(0.4))<1e-9, '');

er = target_state_estimator(est_end{1}, t_sim(1), M(1).r(1), M(1).th(1), true);
checks = add_check(checks, 'Reset when time goes backwards (re-initialises)', ...
    er.n_updates==1 && er.t_init==t_sim(1) && max(abs(er.s - H(1).s(1,:)')) < 1e-12, '');
eq = target_state_estimator();  eq.q_jerk = 77;
eq = target_state_estimator(eq, 1.0, 1500, 0.4, true);
eq = target_state_estimator(eq, 0.5, 1500, 0.4, true);
checks = add_check(checks, 'Reset keeps tuning parameters', ...
    eq.q_jerk == 77 && eq.n_updates == 1 && eq.t_init == 0.5, '');

j_dup = 51;                                        % t = 0.50 s, a radar instant
e1 = target_state_estimator();
for j = 1:j_dup
    e1 = target_state_estimator(e1, t_sim(j), M(2).r(j), M(2).th(j), M(2).valid(j));
end
e2 = target_state_estimator(e1, t_sim(j_dup), M(2).r(j_dup), M(2).th(j_dup), M(2).valid(j_dup));
checks = add_check(checks, 'Repeated call at same t not double-counted', ...
    M(2).valid(j_dup) && isequal(e1.s, e2.s) && e2.n_updates==e1.n_updates && ~e2.updated, '');

ew = target_state_estimator();
ew = target_state_estimator(ew, 0.00, 1000, pi-1e-3, true);
p_w0 = ew.s(1:2);
ew = target_state_estimator(ew, 0.05, 1000, -pi+1e-3, true);   % same place, wrapped angle
checks = add_check(checks, 'Angle wrapping: +-pi crossing does not break the update', ...
    norm(ew.s(1:2) - p_w0) < 50, sprintf('state moved %.2f m', norm(ew.s(1:2)-p_w0)));

idx = 1:5:nS;                                       % same data, 50 ms steps
h05 = run_filter(t_sim(idx), M(2).r(idx), M(2).th(idx), M(2).valid(idx), NaN, false);
dmax = max(max(abs(h05.s - H(2).s(idx,:))));
checks = add_check(checks, 'Step-size invariance: 10 ms vs 50 ms calls (states)', ...
    dmax < 1e-6, sprintf('max state difference %.2e', dmax));
checks = add_check(checks, 'Step-size invariance: update counts equal', ...
    isequal(h05.nupd, H(2).nupd(idx)), '');

pred_u = target_predictor(e1);
e_prop = target_state_estimator(e1, e1.t + P.t_p, NaN, NaN, false);
checks = add_check(checks, 'Predictor state = estimator propagated by Tp', ...
    max(abs(pred_u.s - e_prop.s)) < 1e-9, '');
checks = add_check(checks, 'Predictor covariance = estimator propagated by Tp', ...
    max(max(abs(pred_u.P - e_prop.P))) < 1e-9*max(1,max(abs(e_prop.P(:)))), '');

files = {'target_state_estimator','target_predictor'};
forbid = {'target_truth','target_motion_generator','scnB','scnC'};
for k = 1:2
    txt = fileread(which(files{k}));
    hit = '';
    for q = 1:numel(forbid)
        if contains(txt, forbid{q}), hit = [hit forbid{q} ' ']; end %#ok<AGROW>
    end
    checks = add_check(checks, ['Source of ' files{k} ' has no truth/scenario references'], isempty(hit), hit);
end

npass = sum(cell2mat(checks(:,2)));
fprintf('\n%d / %d hard checks passed.\n', npass, size(checks,1));
if npass == size(checks,1), fprintf('ALL HARD CHECKS PASSED\n'); else, fprintf('SOME HARD CHECKS FAILED - see [FAIL] lines\n'); end

%% 5) PERFORMANCE DIAGNOSTICS (no pass/fail; truth used only here)
rmsf = @(v) sqrt(mean(v.^2));
win  = t_sim >= 1.0;                                          % settled window
pwin = win & (t_sim <= P.t_end - P.t_p + 1e-9);               % prediction window
fprintf('\n---- PERFORMANCE DIAGNOSTICS (not requirements) ----\n');
fprintf('Settled window: t >= 1.0 s. Prediction compared with TRUE position at t + Tp (Tp = %.2f s),\n', P.t_p);
fprintf('so prediction metrics use 1.0 <= t <= %.1f s. R_eng = %.0f m is shown only as a reference.\n\n', P.t_end-P.t_p, P.R_eng);
for s = 1:3
    h = H(s);  e_ = T(s);
    fprintf('Scenario %s\n', scn{s});
    fprintf('  position error     [m]    : RMS %6.2f  max %6.2f   (first second: max %6.2f)\n', ...
        rmsf(e_.pos_err(win)), max(e_.pos_err(win)), max(e_.pos_err(~win)));
    fprintf('  velocity error     [m/s]  : RMS %6.2f  max %6.2f\n', rmsf(e_.vel_err(win)), max(e_.vel_err(win)));
    fprintf('  accel error        [m/s^2]: RMS %6.2f  max %6.2f\n', rmsf(e_.acc_err(win)), max(e_.acc_err(win)));
    fprintf('  predicted position [m]    : RMS %6.2f  max %6.2f\n', rmsf(e_.pred_err(pwin)), max(e_.pred_err(pwin)));
    fprintf('  predicted elevation[deg]  : RMS %6.4f  max %6.4f\n', ...
        rmsf(e_.dth(pwin))*180/pi, max(abs(e_.dth(pwin)))*180/pi);
    fprintf('  miss-equivalent    [m]    : RMS %6.2f  max %6.2f  fraction <= R_eng: %.2f\n', ...
        rmsf(e_.miss(pwin)), max(e_.miss(pwin)), mean(e_.miss(pwin) <= P.R_eng));
    fprintf('  mean NIS (2 dof, ~2 if consistent): %.2f\n\n', mean(h.nis(win), 'omitnan'));
end

wb = t_sim >= 1.0 & t_sim <= 2.0;
fprintf('Blackout effect: B and C share trajectory AND noise; only the blackout differs. Window 1.0-2.0 s.\n');
for s = 2:3
    fprintf('  %s: max predicted-position error %.2f m, mean %.2f m, estimate pos sigma at t=%.2f s: %.2f m\n', ...
        scn{s}, max(T(s).pred_err(wb)), mean(T(s).pred_err(wb)), t_sim(i_last), sqrt(H(s).Ppos(i_last)));
end

q_list = [1 3 10 30 100];
fprintf('\nq sweep (diagnostic): predicted-position RMS [m] / miss-equivalent RMS [m], window 1.0-%.1f s\n', P.t_end-P.t_p);
for qi = 1:numel(q_list)
    fprintf('  q = %5.1f :', q_list(qi));
    for s = 1:3
        hq = run_filter(t_sim, M(s).r, M(s).th, M(s).valid, q_list(qi), false);
        pe = hypot(hq.px - T(s).xf, hq.py - T(s).yf);
        ms = hypot(T(s).xf, T(s).yf) .* abs(sin(wrap_pi(hq.pth - atan2(T(s).yf, T(s).xf))));
        fprintf('   %s: %5.2f / %5.2f', scn{s}, rmsf(pe(pwin)), rmsf(ms(pwin)));
    end
    if q_list(qi) == 10, fprintf('   <- default'); end
    fprintf('\n');
end

%% 6) PLOTS
tp_ok = t_sim <= P.t_end - P.t_p + 1e-9;
tb0 = P.scnC_t_bo0;  tb1 = tb0 + P.scnC_t_bo;

figure('Name','Position: true vs estimated','Position',[30 30 1300 700]);
for s = 1:3
    subplot(2,3,s); hold on; grid on;
    plot(T(s).x, T(s).y, 'k-', 'LineWidth', 1.5);
    plot(H(s).s(1:10:end,1), H(s).s(1:10:end,2), 'r.', 'MarkerSize', 8);
    xlabel('x [m]'); ylabel('y [m]'); title(['Scenario ' scn{s} ': path']);
    legend('true','estimated (every 0.1 s)','Location','best');
    subplot(2,3,3+s); hold on; grid on;
    plot(t_sim, H(s).s(:,1)-T(s).x, 'b');
    plot(t_sim, H(s).s(:,2)-T(s).y, 'r');
    plot(t_sim, T(s).pos_err, 'k', 'LineWidth', 1.2);
    if s == 3, shade_bo(tb0, tb1); end
    xlabel('t [s]'); ylabel('estimate - true [m]'); title(['Scenario ' scn{s} ': position error']);
    legend('x error','y error','|error|','Location','best');
end

figure('Name','Velocity and acceleration estimates','Position',[60 60 1300 700]);
for s = 1:3
    subplot(2,3,s); hold on; grid on;
    plot(t_sim, T(s).vx, 'k-');  plot(t_sim, H(s).s(:,3), 'b-');
    plot(t_sim, T(s).vy, 'k--'); plot(t_sim, H(s).s(:,4), 'r-');
    if s == 3, shade_bo(tb0, tb1); end
    xlabel('t [s]'); ylabel('velocity [m/s]'); title(['Scenario ' scn{s} ': velocity']);
    legend('vx true','vx est','vy true','vy est','Location','best');
    subplot(2,3,3+s); hold on; grid on;
    plot(t_sim, T(s).ax, 'k-');  plot(t_sim, H(s).s(:,5), 'b-');
    plot(t_sim, T(s).ay, 'k--'); plot(t_sim, H(s).s(:,6), 'r-');
    if s == 3, shade_bo(tb0, tb1); end
    xlabel('t [s]'); ylabel('acceleration [m/s^2]'); title(['Scenario ' scn{s} ': acceleration']);
    legend('ax true','ax est','ay true','ay est','Location','best');
end

figure('Name','Prediction error (0.5 s ahead)','Position',[90 90 1300 700]);
for s = 1:3
    subplot(2,3,s); hold on; grid on;
    plot(t_sim(tp_ok), T(s).pred_err(tp_ok), 'b-', 'LineWidth', 1.2);
    plot(t_sim(tp_ok), T(s).miss(tp_ok), 'r-');
    plot([0 P.t_end-P.t_p], [P.R_eng P.R_eng], 'k:');
    if s == 3, shade_bo(tb0, tb1); end
    xlabel('t [s]'); ylabel('[m]'); title(['Scenario ' scn{s} ': prediction error']);
    legend('|pred - true pos at t+Tp|','miss-equivalent','R_{eng} (reference)','Location','best');
    subplot(2,3,3+s); hold on; grid on;
    plot(t_sim(tp_ok), T(s).dth(tp_ok)*180/pi, 'b-');
    plot(t_sim(tp_ok),  H(s).psth(tp_ok)*180/pi, 'k--');
    plot(t_sim(tp_ok), -H(s).psth(tp_ok)*180/pi, 'k--');
    if s == 3, shade_bo(tb0, tb1); end
    xlabel('t [s]'); ylabel('[deg]'); title(['Scenario ' scn{s} ': predicted elevation error']);
    legend('error','\pm1\sigma from filter','Location','best');
end

figure('Name','Scenario C blackout vs B','Position',[120 120 1100 650]);
w = t_sim >= 0.5 & t_sim <= 3.0 & tp_ok;
subplot(2,1,1); hold on; grid on;
plot(t_sim(w), T(2).pred_err(w), 'b-', 'LineWidth', 1.3);
plot(t_sim(w), T(3).pred_err(w), 'r-', 'LineWidth', 1.3);
shade_bo(tb0, tb1);
xlabel('t [s]'); ylabel('prediction position error [m]');
title('B (no blackout) vs C (blackout 1.0-1.5 s): same trajectory and noise');
legend('B','C','Location','best');
subplot(2,1,2); hold on; grid on;
plot(t_sim(w), sqrt(H(2).Ppos(w)), 'b-', 'LineWidth', 1.3);
plot(t_sim(w), sqrt(H(3).Ppos(w)), 'r-', 'LineWidth', 1.3);
shade_bo(tb0, tb1);
xlabel('t [s]'); ylabel('estimate position sigma [m]');
title('Filter-reported position uncertainty sqrt(trace P_{pos})');
legend('B','C','Location','best');

%% Local functions (must stay at the end of the script)
function [h, est] = run_filter(t, r, th, valid, q_override, full)
    n = numel(t);
    est = target_state_estimator();
    if ~isnan(q_override), est.q_jerk = q_override; end
    h.init = false(n,1);  h.s = nan(n,6);  h.Ppos = nan(n,1);  h.Pfin = false(n,1);
    h.P_ok = false(n,1);  h.t_est = nan(n,1);  h.nupd = nan(n,1);  h.updated = false(n,1);
    h.nis = nan(n,1);  h.tinit = nan(n,1);
    h.px = nan(n,1);  h.py = nan(n,1);  h.pth = nan(n,1);  h.pvalid = false(n,1);
    h.psth = nan(n,1);  h.ptt = nan(n,1);
    for j = 1:n
        est  = target_state_estimator(est, t(j), r(j), th(j), valid(j));
        pred = target_predictor(est);
        h.init(j)    = est.initialized;
        h.s(j,:)     = est.s';
        h.t_est(j)   = est.t;
        h.nupd(j)    = est.n_updates;
        h.updated(j) = est.updated;
        h.tinit(j)   = est.t_init;
        h.Pfin(j)    = all(isfinite(est.P(:)));
        if h.Pfin(j)
            h.Ppos(j) = trace(est.P(1:2,1:2));
            if full
                Pj = est.P;
                ev = eig(0.5*(Pj+Pj'));
                h.P_ok(j) = norm(Pj-Pj','fro') <= 1e-9*max(1,norm(Pj,'fro')) && ...
                            min(ev) >= -1e-9*max(abs(ev));
            else
                h.P_ok(j) = true;
            end
        end
        if est.updated && est.n_updates > 1, h.nis(j) = est.nis; end
        h.pvalid(j) = pred.valid;
        h.px(j) = pred.x_p;  h.py(j) = pred.y_p;  h.pth(j) = pred.theta_p;
        h.psth(j) = pred.sigma_theta;  h.ptt(j) = pred.t_target;
    end
end

function checks = add_check(checks, name, ok, detail)
    ok = all(logical(ok(:)));
    if ok, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('[%s] %-66s %s\n', tag, name, detail);
    checks(end+1,:) = {name, ok, detail};
end

function a = wrap_pi(a)
    a = atan2(sin(a), cos(a));
end

function shade_bo(t0, t1)
    yl = ylim;
    hp = patch([t0 t1 t1 t0], [yl(1) yl(1) yl(2) yl(2)], [1 0.85 0.85], ...
               'EdgeColor','none', 'FaceAlpha',0.5, 'HandleVisibility','off');
    uistack(hp, 'bottom');
end