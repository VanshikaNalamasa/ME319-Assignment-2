% test_target_motion_generator.m
% Independent checks + plots for target_truth.m and target_motion_generator.m
% Requires params.m, target_truth.m, target_motion_generator.m in the same folder.
clear; clc; close all;

P   = params();
Ts  = P.Ts_radar;
scn = {'A','B','C'};
checks = cell(0,3);
fprintf('=== test_target_motion_generator ===\n\n');

%% 1) Dense ground truth, 1 kHz, 0..10 s
t_f = (0:0.001:P.t_end)';
for s = 1:3
    [x,y,vx,vy,ax,ay] = target_truth(t_f, scn{s});
    D(s).x = x;  D(s).y = y;  D(s).vx = vx;  D(s).vy = vy;  D(s).ax = ax;  D(s).ay = ay;
    D(s).r = hypot(x,y);  D(s).th = atan2(y,x);
    D(s).speed = hypot(vx,vy);  D(s).acc = hypot(ax,ay);
end
checks = add_check(checks, 'Time grid covers 0..10 s', ...
    numel(t_f)==10001 && abs(t_f(end)-P.t_end)<1e-9, sprintf('%d points', numel(t_f)));

%% 2) Target limits
for s = 1:3
    checks = add_check(checks, ['Range within [800,2000] m, scenario ' scn{s}], ...
        min(D(s).r) >= P.r_min && max(D(s).r) <= P.r_max, ...
        sprintf('%.1f ... %.1f m', min(D(s).r), max(D(s).r)));
    checks = add_check(checks, ['Speed <= 150 m/s, scenario ' scn{s}], ...
        max(D(s).speed) <= P.v_max, ...
        sprintf('%.2f ... %.2f m/s', min(D(s).speed), max(D(s).speed)));
    checks = add_check(checks, ['Accel <= 20 m/s^2, scenario ' scn{s}], ...
        max(D(s).acc) <= P.a_max, sprintf('max %.3f m/s^2', max(D(s).acc)));
    checks = add_check(checks, ['Elevation within turret limits, scenario ' scn{s}], ...
        min(D(s).th) >= P.theta_min_rad && max(D(s).th) <= P.theta_max_rad, ...
        sprintf('%.2f ... %.2f deg', min(D(s).th)*180/pi, max(D(s).th)*180/pi));
end

%% 3) Kinematic consistency (v = dx/dt, a = dv/dt) and smoothness
tt = (0.01:0.01:9.99)';  h = 1e-3;
for s = 1:3
    [xp,yp,vxp,vyp] = target_truth(tt+h, scn{s});
    [xm,ym,vxm,vym] = target_truth(tt-h, scn{s});
    [~,~,vx0,vy0,ax0,ay0] = target_truth(tt, scn{s});
    err_v = max(abs([(xp-xm)/(2*h) - vx0; (yp-ym)/(2*h) - vy0]));
    err_a = max(abs([(vxp-vxm)/(2*h) - ax0; (vyp-vym)/(2*h) - ay0]));
    checks = add_check(checks, ['v consistent with d(pos)/dt, scenario ' scn{s}], ...
        err_v < 1e-4, sprintf('max error %.2e m/s', err_v));
    checks = add_check(checks, ['a consistent with d(vel)/dt, scenario ' scn{s}], ...
        err_a < 1e-4, sprintf('max error %.2e m/s^2', err_a));

    dpos = max(hypot(diff(D(s).x), diff(D(s).y)));
    dvel = max(hypot(diff(D(s).vx), diff(D(s).vy)));
    dacc = max(hypot(diff(D(s).ax), diff(D(s).ay)));
    checks = add_check(checks, ['No jumps in pos/vel/acc per 1 ms, scenario ' scn{s}], ...
        dpos < 0.2 && dvel < 0.1 && dacc < 0.1, ...
        sprintf('max step: %.3f m, %.4f m/s, %.4f m/s^2', dpos, dvel, dacc));
end

%% 4) Scenario structure and maneuver onset
pre = t_f <= P.scnB_maneuver_start;
checks = add_check(checks, 'Scenario A has zero acceleration', ...
    max(D(1).acc) < 1e-12, sprintf('max %.2e', max(D(1).acc)));
for s = 2:3
    a_pre  = max(D(s).acc(pre));
    idx    = find(D(s).acc > 0.01 & t_f > P.scnB_maneuver_start, 1, 'first');
    t_on   = t_f(idx);
    a_post = max(D(s).acc(t_f > P.scnB_maneuver_start));
    checks = add_check(checks, ['Zero acceleration before t=1 s, scenario ' scn{s}], ...
        a_pre < 1e-12, sprintf('max %.2e', a_pre));
    checks = add_check(checks, ['Maneuver onset at t ~ 1 s, scenario ' scn{s}], ...
        ~isempty(idx) && t_on >= 1.0 && t_on <= 1.05, ...
        sprintf('|a| first exceeds 0.01 m/s^2 at t = %.3f s', t_on));
    checks = add_check(checks, ['Real maneuver after onset, scenario ' scn{s}], ...
        a_post > 5, sprintf('peak %.2f m/s^2', a_post));
    dd = max(abs([D(s).x(pre)-D(1).x(pre); D(s).y(pre)-D(1).y(pre); ...
                  D(s).vx(pre)-D(1).vx(pre); D(s).vy(pre)-D(1).vy(pre)]));
    checks = add_check(checks, ['Identical to A before t=1 s, scenario ' scn{s}], ...
        dd < 1e-9, sprintf('max diff %.2e', dd));
end
checks = add_check(checks, 'Scenario C trajectory equals B (assumption)', ...
    max(abs([D(3).x-D(2).x; D(3).y-D(2).y])) < 1e-12, '');

%% 5) Radar sampling at 20 Hz via the generator
k  = 0:200;  ts = k*Ts;
R  = nan(3,numel(ts));  TH = R;  V = false(3,numel(ts));
for s = 1:3
    for j = 1:numel(ts)
        [R(s,j), TH(s,j), V(s,j)] = target_motion_generator(ts(j), scn{s});
    end
end
bo = (ts >= P.scnC_t_bo0 - 1e-9) & (ts < P.scnC_t_bo0 + P.scnC_t_bo - 1e-9);

checks = add_check(checks, 'A: all 201 samples valid', all(V(1,:)), sprintf('%d valid', nnz(V(1,:))));
checks = add_check(checks, 'B: all 201 samples valid', all(V(2,:)), sprintf('%d valid', nnz(V(2,:))));
checks = add_check(checks, 'C: valid flags match blackout [1.0,1.5) s', ...
    isequal(V(3,:), ~bo), sprintf('%d valid, %d invalid', nnz(V(3,:)), nnz(~V(3,:))));
idx_inv = find(~V(3,:));
checks = add_check(checks, 'C: blackout starts at 1.00 s', abs(ts(idx_inv(1))-1.00)<1e-9, ...
    sprintf('first invalid t = %.2f s', ts(idx_inv(1))));
checks = add_check(checks, 'C: last blackout sample 1.45 s, 10 samples, contiguous', ...
    abs(ts(idx_inv(end))-1.45)<1e-9 && numel(idx_inv)==10 && all(diff(idx_inv)==1), ...
    sprintf('last invalid t = %.2f s', ts(idx_inv(end))));
checks = add_check(checks, 'C: measurements resume at 1.50 s', ...
    V(3, idx_inv(end)+1) && abs(ts(idx_inv(end)+1)-1.50)<1e-9, ...
    sprintf('resume t = %.2f s', ts(idx_inv(end)+1)));
for s = 1:3
    checks = add_check(checks, ['NaN exactly when invalid, scenario ' scn{s}], ...
        all(isnan(R(s,~V(s,:)))) && all(isnan(TH(s,~V(s,:)))) && ...
        all(isfinite(R(s,V(s,:)))) && all(isfinite(TH(s,V(s,:)))), '');
end

%% 6) Validity on the dense 1 ms grid (valid only on 50 ms instants)
Vd = false(3,numel(t_f));
for s = 1:3
    for j = 1:numel(t_f)
        [~,~,Vd(s,j)] = target_motion_generator(t_f(j), scn{s});
    end
end
checks = add_check(checks, 'A: valid exactly at 201 radar instants (1 ms calls)', ...
    nnz(Vd(1,:))==201 && max(abs(t_f(Vd(1,:)) - ts')) < 1e-9, sprintf('%d valid', nnz(Vd(1,:))));
checks = add_check(checks, 'C: valid exactly at 191 radar instants (1 ms calls)', ...
    nnz(Vd(3,:))==191 && max(abs(t_f(Vd(3,:)) - ts(~bo)')) < 1e-9, sprintf('%d valid', nnz(Vd(3,:))));
[r1,th1,v1] = target_motion_generator(0.0251, 'A');
checks = add_check(checks, 'Off-instant call (t=0.0251) gives no measurement', ...
    ~v1 && isnan(r1) && isnan(th1), '');
[~,~,v2] = target_motion_generator(0.05 + 1e-9, 'A');
[~,~,v3] = target_motion_generator(0.05 + 1e-3, 'A');
checks = add_check(checks, 'Tolerance: 1 ns off valid, 1 ms off invalid', v2 && ~v3, '');

%% 7) Measurements vs truth, reproducibility, interface
RT = nan(3,numel(ts));  THT = RT;
for s = 1:3
    [xs,ys] = target_truth(ts', scn{s});
    RT(s,:) = hypot(xs,ys)';  THT(s,:) = atan2(ys,xs)';
end
RES_r = R - RT;  RES_t = TH - THT;
for s = 1:3
    vv = V(s,:);
    checks = add_check(checks, ['Residuals within 6 sigma, scenario ' scn{s}], ...
        max(abs(RES_r(s,vv))) < 6*P.sigma_r && max(abs(RES_t(s,vv))) < 6*P.sigma_theta_rad, ...
        sprintf('max |dr| = %.2f m, max |dtheta| = %.2e rad', max(abs(RES_r(s,vv))), max(abs(RES_t(s,vv)))));
end
both = V(1,:) & V(2,:);
checks = add_check(checks, 'Same noise realisation in A and B', ...
    max(abs(RES_r(1,both)-RES_r(2,both))) < 1e-9 && max(abs(RES_t(1,both)-RES_t(2,both))) < 1e-12, '');
both = V(2,:) & V(3,:);
checks = add_check(checks, 'B and C measurements identical outside blackout', ...
    max(abs(R(2,both)-R(3,both))) < 1e-12 && max(abs(TH(2,both)-TH(3,both))) < 1e-15, '');

a1 = cell(1,3);  a2 = cell(1,3);
[a1{:}] = target_motion_generator(1.2, 'B');
[a2{:}] = target_motion_generator(1.2, 'B');
checks = add_check(checks, 'Repeated call at same t gives identical output', isequaln(a1,a2), '');

Rrev = nan(1,numel(ts));
for j = numel(ts):-1:1
    Rrev(j) = target_motion_generator(ts(j), 'C');
end
checks = add_check(checks, 'Output independent of call order (reverse sweep, C)', isequaln(Rrev, R(3,:)), '');

rng(123);  u1 = rand(1,5);
rng(123);  target_motion_generator(0.5,'B');  target_motion_generator(0.55,'B');  u2 = rand(1,5);
checks = add_check(checks, 'Global random-number stream not consumed', isequal(u1,u2), '');

[ra,tha,va] = target_motion_generator(0.5,'B');
[rb,thb,vb] = target_motion_generator(0.5, 2);
checks = add_check(checks, 'Scenario given as char or number gives same result', ...
    ra==rb && tha==thb && va==vb, '');
[r0,th0,v0] = target_motion_generator(0.5,'A');
checks = add_check(checks, 'Output sizes/types (scalar, logical flag)', ...
    isscalar(r0) && isscalar(th0) && isscalar(v0) && islogical(v0), '');
bad = false;
try, target_motion_generator(0,'Z'); catch, bad = true; end
checks = add_check(checks, 'Unknown scenario raises an error', bad, '');

%% 8) Noise statistics over 4000 repeated measurements (Scenario A)
N  = 4000;  tn = (0:N-1)'*Ts;
rn = nan(N,1);  thn = nan(N,1);  vn = false(N,1);
for j = 1:N
    [rn(j), thn(j), vn(j)] = target_motion_generator(tn(j), 'A');
end
[xn,yn] = target_truth(tn, 'A');
er = rn - hypot(xn,yn);
et = thn - atan2(yn,xn);
er0 = er - mean(er);  et0 = et - mean(et);
rho  = sum(er0.*et0) / sqrt(sum(er0.^2)*sum(et0.^2));
lag_r = sum(er0(1:end-1).*er0(2:end)) / sum(er0.^2);
lag_t = sum(et0(1:end-1).*et0(2:end)) / sum(et0.^2);
lim = 4/sqrt(N);

checks = add_check(checks, 'Noise: all 4000 samples valid', all(vn), '');
checks = add_check(checks, 'Range noise mean ~ 0', abs(mean(er)) < 4*P.sigma_r/sqrt(N), ...
    sprintf('mean = %+.4f m', mean(er)));
checks = add_check(checks, 'Range noise std ~ 1 m', abs(std(er)/P.sigma_r - 1) < 0.05, ...
    sprintf('std = %.4f m', std(er)));
checks = add_check(checks, 'Angle noise mean ~ 0', abs(mean(et)) < 4*P.sigma_theta_rad/sqrt(N), ...
    sprintf('mean = %+.2e rad', mean(et)));
checks = add_check(checks, 'Angle noise std ~ 0.02 deg (units correct)', ...
    abs(std(et)/P.sigma_theta_rad - 1) < 0.05, ...
    sprintf('std = %.5f deg = %.3e rad', std(et)*180/pi, std(et)));
checks = add_check(checks, 'Range/angle noise uncorrelated', abs(rho) < lim, sprintf('rho = %+.4f', rho));
checks = add_check(checks, 'No sample-to-sample correlation (lag 1)', ...
    abs(lag_r) < lim && abs(lag_t) < lim, sprintf('lag-1: %+.4f (r), %+.4f (theta)', lag_r, lag_t));

%% Summary
npass = sum(cell2mat(checks(:,2)));
fprintf('\n%d / %d checks passed.\n', npass, size(checks,1));
if npass == size(checks,1), fprintf('ALL CHECKS PASSED\n'); else, fprintf('SOME CHECKS FAILED - see [FAIL] lines\n'); end

%% Plots
col = {'b','r','g'};  lstyle = {'-','-','--'};
figure('Name','Target truth','Position',[40 40 1250 700]);

subplot(2,3,1); hold on; grid on;
h = gobjects(1,3);
for s = 1:3, h(s) = plot(D(s).x, D(s).y, [col{s} lstyle{s}], 'LineWidth', 1.5); end
phi = linspace(0, pi/2, 100);
plot(P.r_min*cos(phi), P.r_min*sin(phi), 'k:', 'HandleVisibility','off');
plot(P.r_max*cos(phi), P.r_max*sin(phi), 'k:', 'HandleVisibility','off');
plot(0, 0, 'k^', 'MarkerFaceColor','k', 'HandleVisibility','off');
axis equal; xlabel('x [m]'); ylabel('y [m]'); title('Trajectory (dotted: 800 m, 2000 m)');
legend(h, {'A','B','C (= B)'}, 'Location','best');

subplot(2,3,2); hold on; grid on;
for s = 1:3, plot(t_f, D(s).r, [col{s} lstyle{s}], 'LineWidth', 1.5); end
plot([0 P.t_end],[P.r_min P.r_min],'k:'); plot([0 P.t_end],[P.r_max P.r_max],'k:');
xlabel('t [s]'); ylabel('range [m]'); title('Target range');

subplot(2,3,3); hold on; grid on;
for s = 1:3, plot(t_f, D(s).th*180/pi, [col{s} lstyle{s}], 'LineWidth', 1.5); end
xlabel('t [s]'); ylabel('elevation [deg]'); title('Target elevation (turret limits -10..85)');

subplot(2,3,4); hold on; grid on;
for s = 1:3, plot(t_f, D(s).speed, [col{s} lstyle{s}], 'LineWidth', 1.5); end
plot([0 P.t_end],[P.v_max P.v_max],'k:');
xlabel('t [s]'); ylabel('speed [m/s]'); title('Target speed (limit 150)');

subplot(2,3,5); hold on; grid on;
for s = 1:3, plot(t_f, D(s).acc, [col{s} lstyle{s}], 'LineWidth', 1.5); end
plot([0 P.t_end],[P.a_max P.a_max],'k:');
xlabel('t [s]'); ylabel('|a| [m/s^2]'); title('Target acceleration (limit 20)');

subplot(2,3,6); hold on; grid on;
for s = 1:3, stairs(ts, double(V(s,:)) + 2*(3-s), col{s}, 'LineWidth', 1.2); end
set(gca,'YTick',0:5,'YTickLabel',{'C:0','C:1','B:0','B:1','A:0','A:1'});
ylim([-0.5 5.5]); xlabel('t [s]'); title('meas\_valid at 20 Hz');
saveas(gcf, 'tmg_fig1_truth.png');

figure('Name','Scenario C measurements','Position',[80 80 1250 700]);
w_f = t_f >= 0.5 & t_f <= 2.5;   w_s = ts >= 0.5 & ts <= 2.5;
ok_s = w_s & V(3,:);             bo_s = w_s & ~V(3,:);

subplot(2,2,1); hold on; grid on;
plot(t_f(w_f), D(3).r(w_f), 'k-');
plot(ts(ok_s), R(3,ok_s), 'b.', 'MarkerSize', 10);
plot(ts(bo_s), RT(3,bo_s), 'rx');
xlim([0.5 2.5]); shade_bo(P);
xlabel('t [s]'); ylabel('range [m]'); title('Scenario C: range (shaded = blackout)');
legend('true','measured','no measurement','Location','best');

subplot(2,2,2); hold on; grid on;
plot(t_f(w_f), D(3).th(w_f)*180/pi, 'k-');
plot(ts(ok_s), TH(3,ok_s)*180/pi, 'b.', 'MarkerSize', 10);
plot(ts(bo_s), THT(3,bo_s)*180/pi, 'rx');
xlim([0.5 2.5]); shade_bo(P);
xlabel('t [s]'); ylabel('elevation [deg]'); title('Scenario C: elevation');

subplot(2,2,3); histogram(er, 40); grid on;
xlabel('range error [m]'); title(sprintf('Range noise, std = %.3f m', std(er)));
subplot(2,2,4); histogram(et*180/pi, 40); grid on;
xlabel('elevation error [deg]'); title(sprintf('Angle noise, std = %.5f deg', std(et)*180/pi));
saveas(gcf, 'tmg_fig2_measurements.png');

%% Local functions (must stay at the end of the script)
function checks = add_check(checks, name, ok, detail)
    ok = all(logical(ok(:)));
    if ok, tag = 'PASS'; else, tag = 'FAIL'; end
    fprintf('[%s] %-62s %s\n', tag, name, detail);
    checks(end+1,:) = {name, ok, detail};
end

function shade_bo(P)
    yl = ylim;
    t0 = P.scnC_t_bo0;  t1 = t0 + P.scnC_t_bo;
    hp = patch([t0 t1 t1 t0], [yl(1) yl(1) yl(2) yl(2)], [1 0.85 0.85], ...
               'EdgeColor','none', 'FaceAlpha',0.5, 'HandleVisibility','off');
    uistack(hp, 'bottom');
end