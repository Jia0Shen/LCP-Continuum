clear; clc

ctr = build_ctr();

ctr_est = ctr;
ctr_est.E = ctr.E * (1 + 0.0);  % 5% uncertainty
% ctr_est.G = ctr.G * (1 + 0);  % 5% uncertainty
% ctr_est.E = ctr.E * (1 + (rand(1)-0.5)/10);  % 5% uncertainty
% ctr_est.G = ctr.G * (1 + (rand(1)-0.5)/10);  % 5% uncertainty
ctr_est.K1 = diag([ctr_est.E*ctr.I1, ctr_est.E*ctr.I1, ctr_est.G*ctr.J1]);
ctr_est.K2 = diag([ctr_est.E*ctr.I2, ctr_est.E*ctr.I2, ctr_est.G*ctr.J2]);
ctr_est.K3 = diag([ctr_est.E*ctr.I3, ctr_est.E*ctr.I3, ctr_est.G*ctr.J3]);
% change precurvature
% ctr.ka1 = 5; ctr.ka2 = 7; ctr.ka3 = 4;
%% Create Obstacles, all in meter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% plane_fl = createPlane(p0, n, 1e-5);

q0 = [-3/2*pi -0.1 -3/4*pi 0 0 0]';
q1 = [-3/2*pi -0.0 -3/4*pi 0 0 0]';
% q2 = [2*pi 0.02 pi 0 2*pi 0]';
dalpha = deg2rad(1);   % rotate in 1 degree per step
dbeta = 2e-5;   % move in 2e-2 mm per step

q_traj1 = LinearInterplNdim([q0,q1], dbeta);

w0 = [0 0 0 0 0 0]';
yu0_guess = zeros(5,1);

%% 
n1 = size(q_traj1, 2);

%start 
tic

tipContact = false;
% yu0 = yu0_guess;

tip_traj_ini = [];

state_traj_ini = [];

[~, p_ini, yu0, b_res, g_ini] = three_tube_fk(ctr, q0, w0, yu0_guess);
[~, p1, ~, b_res1, g_ini1] = three_tube_fk(ctr, q1, w0, yu0_guess);


%  def plane
wall_dis = p_ini{3}(end,3)+0.001;
p0 = [0;0;wall_dis];
n = [0;0;-1];
mu = 0.2;

plane = createPlane(p0, n, mu);

plane_fl = plane;
plane_fl.mu = 0;

% visualize

% figure()
% hold on
% plotConfig3D([640 640])
% xlabel('x (m)')
% ylabel('y (m)')
% zlabel('z (m)')
% axis equal
% axislim = 0.15;
% xlim([-axislim, axislim])
% ylim([-axislim, axislim])
% % xlim([-axislim, 0.02])
% % ylim([-0.04, axislim])
% zlim([-0.02, 0.52])
% 
% view(30,45)
% h = plotPlane(plane, 0.5);
% h.FaceColor = 'k';
% h.FaceAlpha = 0.2;
% simple_plot_ctr(g_ini, 1, 'r');
% simple_plot_ctr(g_ini1, 1, 'r');


% initialize
state_ini.p = p_ini;
state_ini.w = w0;
state_ini.yu0 = yu0;
state_ini.g = g_ini;
state_ini.q = q0;
state_prev = state_ini;

d = 60;


% start the trajecotry
% save 
loadStates = true;
if loadStates
    load("simu_out\ctr_pdi_ini.mat");
    state_st = state_traj_ini(end);
else

    for i = 1:n1
    
        qi = q_traj1(:, i);
    
        if tipContact
            % assume obs not moving here
            [state_update,tipContact] = getFrictionalTipContactCTR(ctr, plane_fl, qi, state_prev, d);
    
            state_prev = state_update;
    
            p_tip = state_update.p{3}(end,:)';
    
        else
            [~, p, yu0, b_res, g] = three_tube_fk(ctr, qi, w0, yu0);
    
            state_update.p = p;
            state_update.w = w0;
            state_update.yu0 = yu0;
            state_update.q = qi;
            state_update.g = g;
    
            state_prev = state_update;
    
            % check if tip contact
            p_tip = p{3}(end,:)';
    
            if plane.computeDepth(p_tip) >= 0
                tipContact = true;
            else
                tipContact = false;
            end
    
        end
    
        tip_traj_ini = [tip_traj_ini, p_tip];
        state_traj_ini = [state_traj_ini, state_update];
    
    end
    
    tolTime = toc;
    
    % disp(n/tolTime)
    fprintf('Finished pushing in %.4f s, get normal force. \n', tolTime)

    save("simu_out\ctr_pdi_ini.mat", "state_traj_ini", "tip_traj_ini");
    state_st = state_traj_ini(end);
end

%%

% design path
rp = norm(tip_traj_ini(1:2, end));
x_tip = tip_traj_ini(1, end);
y_tip = tip_traj_ini(2, end);
z_tip = tip_traj_ini(3, end);
theta0 = atan2(y_tip, x_tip);

directSign = 1;   % +1 for CCW, -1 for CW.
% thetas = linspace(theta0, theta0+2*directSign*pi, n_path);
theta1 = theta0+directSign*pi*1/8;
% thetas = LinearInterpl([theta0, theta1, theta0, theta1], 0.002);
thetas = LinearInterpl([theta0, theta1, theta0], 0.0005);
n_path = length(thetas);
x_path = [rp*cos(thetas); rp*sin(thetas); z_tip*ones(1,n_path)];

% initialize
q_st = state_st.q;
w_st = state_st.w;
yu0_st = state_st.yu0;

xc = tip_traj_ini(:,end);
[J, C, ~, ~, B] = cal_Jacobian_ivp_one(ctr, q_st, w_st, yu0_st);
Jc = J(1:3,:);

% gain
KP = 15*diag([1 1 1]);
KI = 20*diag([1 1 1]);
Ks = 0.9;
err_int = zeros(3, 1);
s = 0;

dt = 0.01;  % 100 Hz
tipContact = true;
d = 60;
qc = q_st;
yu0_c = yu0_st;
state_c = state_st;
state_est = state_st;

state_traj = [];
s_list = [];

% change the COF
plane.mu = 0.37;

dq_idx = [0 0 0]';

for i = 1:n_path

    % find jacobian - dq
    xd = x_path(:, i);

    % opt 1: this assumes known tip position control.
    % dx = (xd-xc)/dt;   

    % opt 2: this assumes open-loop control.
    if i == 1
        dx = [0;0;0];
    else 
        dx = (x_path(:, i) - x_path(:, i-1))/dt; 
    end    

    % print 
    if mod(i, 50) == 0; fprintf("Finished %d/%d iterations. \n", i, n_path); end

    [J, C, ~, ~, B] = cal_Jacobian_ivp_one(ctr_est, state_est.q, state_est.w, state_est.yu0);

    % get updated Jacovian using LCP-model

    % J_lcp = [];

    % [J, C, ~, ~, B] = cal_Jacobian_ivp_one(ctr, state_c.q, w_st, yu0_c);

    % use all 6 DOF for j_idx = 1:6;
    j_idx = [1 3 5];
    % j_idx = 1:6;
    Jv = J(1:3,j_idx);

    Jv_inv = Jv'*(inv(Jv*Jv' + 0.00001 * eye(3)));   % singularity-robust

    % opt1 PID control
    % err_int = err_int + (xd-xc)*dt;
    % dq_idx = Jv_inv*(dx+KP*(xd-xc)+KI*err_int);

    % opt2 Sliding mode control
    sat = @(z, thres) min(max(z, -thres), thres);

    s = Jv * dq_idx - dx + KP*(xc-xd);
    % record s
    s_list = [s_list, s];

    dq_idx = Jv_inv*(dx+KP*(xd-xc)+Ks*sat(s, 0.1));
    
    dq = zeros(6,1);
    dq(j_idx) = dq_idx;

    % dq = zeros(size(qc));
    % dq = [1 0 1 0 1 0]' * 5e-1;

    q_prev = qc;
    qc = q_prev + dq*dt;

    % qc = q_traj2;   % test a round traj

    % estimate shape using (inaccurate) model.
    % [~, p_est, yu0_c, b_est, g_est] = three_tube_fk(ctr, qc, w_st, yu0_c);
    % [~, p_est, yu0_c, b_est, g_est] = three_tube_fk(ctr, qc, state_c.w, state_c.yu0);

    % feed forward simulator (ground truth), noise included?

    state_prev = state_c;
    state_est_prev = state_est;

    % update qc
    % if i >= 198 && i <= 237
    %     q_prev = qc;
    %     qc = stepLCPControl(ctr, plane, xd, state_prev, d);
    %     disp(i)
    % else
    %     q_prev = qc;
    %     qc = q_prev + dq*dt;
    % end



    if tipContact
        % assume obs not moving here
        [state_c,tipContact] = getFrictionalTipContactCTR(ctr, plane, qc, state_prev, d);

        % estimate the state using frictionless model
        [state_est,~] = getFrictionalTipContactCTR(ctr_est, plane_fl, qc, state_est_prev, d);

        % state_prev = state_update;
        xc = state_c.p{3}(end,:)';

    else
        % tip wrench is 0 if no contact.
        [~, p, yu0, b_res, g] = three_tube_fk(ctr, qc, w0, yu0);

        state_c.p = p;
        state_c.w = w0;
        state_c.yu0 = yu0;
        state_c.q = qi;
        state_c.g = g;

        % state_prev = state_update;

        % check if tip contact
        xc = p{3}(end,:)';

        if plane.computeDepth(xc) >= 0
            tipContact = true;
        else
            tipContact = false;
        end

    end

    % update sliding mode control
    % s = Jv * dq_idx - dx + KP*(xc-xd);

    % further update other states.
    % xc = p_tip;
    yu0_c = state_c.yu0;
    state_record = state_c;
    state_record.target = xd;
    state_record.error = norm(xd-xc);
    state_record.g_est = state_est.g;
    state_record.dt = dt;
    state_traj = [state_traj, state_record];

end

%%

tipContact = true;
qc = q_st;
yu0_c = yu0_st;
state_c = state_st;
state_est = state_st;
state_traj_lcp = [];

mpc_horizon = 5;

for i = 1:n_path

    % find jacobian - dq
    xd = x_path(:, i);
    
    % feed forward simulator (ground truth), noise included?

    state_prev = state_c;

    % find qc given the desired next step
    % [state_c,tipContact] = getFrictionalTipContactCTR(ctr, plane, qc, state_prev, d);

    % if i >= 0
    %     qc = stepLCPControl(ctr, plane, xd, state_prev, d);
    % else
    %     qc = state_traj(i).q;
    % end

    qc = stepLCPControl(ctr, plane, xd, state_prev, d);

    % find qc using mpc
    % if i <= n_path - mpc_horizon
    %     xd_mpc = x_path(:, i-1+(1:mpc_horizon));   % i, i+1, ... steps
    %     qc = stepLCPControl_mpc(ctr, plane, xd_mpc, state_prev, d);
    % else
    %     xd_mpc = x_path(:, i:n_path);   % i, i+1, ..., n_path
    %     qc = stepLCPControl_mpc(ctr, plane, xd_mpc, state_prev, d);
    % end

    % compare with frictionless.
    % qc_compare = state_traj(i).q;

    if tipContact
        % assume obs not moving here
        [state_c,tipContact] = getFrictionalTipContactCTR(ctr, plane, qc, state_prev, d);

        % state_prev = state_update;
        xc = state_c.p{3}(end,:)';

    else
        % tip wrench is 0 if no contact.
        [~, p, yu0, b_res, g] = three_tube_fk(ctr, qc, w0, yu0);

        state_c.p = p;
        state_c.w = w0;
        state_c.yu0 = yu0;
        state_c.q = qc;
        state_c.g = g;

        % state_prev = state_update;

        % check if tip contact
        xc = p{3}(end,:)';

        if plane.computeDepth(xc) >= 0
            tipContact = true;
        else
            tipContact = false;
        end

    end

    % further update other states.
    % xc = p_tip;
    yu0_c = state_c.yu0;
    state_record2 = state_c;
    state_record2.target = xd;
    state_record2.error = norm(xd-xc);
    state_record2.g_est = [];
    state_record2.dt = dt;
    state_traj_lcp = [state_traj_lcp, state_record2];

    % print 
    if mod(i, 1) == 0; fprintf("Finished %d/%d iterations, error is %.6f. \n", i, n_path, norm(xd-xc)); end

end

save("simu_out\ctr_pid_traj.mat", 'state_traj_lcp', 'state_traj', 'x_path')

%% plot the traj

% plot_traj = [state_traj_ini, state_traj];
% plot_traj = state_traj_ini;

% down sampled for visualization.
idx_sampled = [1:20:n_path, n_path];

plot_traj = state_traj(idx_sampled);
plot_traj2 = state_traj_lcp(idx_sampled);

% plot_traj = state_traj;

plot_tip_traj = [];
plot_tip_traj2 = [];
n_traj = length(plot_traj);

% close(vid1)

totalTime = 10;

vid1 = VideoWriter('simu_out/vid_ctr_pid2', 'MPEG-4');
vid1.FrameRate = min(round(n_traj/totalTime), 110);
open(vid1);

figure()
hold on
plotConfig3D([640 640])

% plot3Mat(tip_traj)

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis equal
% axislim = 0.2;
% xlim([-axislim, axislim])
% ylim([-axislim, axislim])
xlim([-0.23, 0.02])
ylim([-0.05, 0.2])
zlim([-0.02, 0.5])

view(5,60)
% view(0, 90)

h = plotPlane(plane, 0.5);
h.FaceColor = 'k';
h.FaceAlpha = 0.2;
hs = []; ht = []; hs2 = []; ht2 = [];

% plot path
plot3(x_path(1,idx_sampled), x_path(2,idx_sampled), x_path(3,idx_sampled), 'g--', 'LineWidth',2)


for i = 1:n_traj

    % statei = plot_traj(i);
    statei = plot_traj(i);
    statei2 = plot_traj2(i);

    plot_tip_traj = [plot_tip_traj, statei.p{3}(end,:)'];
    plot_tip_traj2 = [plot_tip_traj2, statei2.p{3}(end,:)'];
    % plot_tip_traj_est = [plot_tip_traj_est, statei.g_est{3}(end, 13:15)'];

    delete([hs hs2 ht ht2])

    hs = simple_plot_ctr(statei.g, 1, 'b');
    hs2 = simple_plot_ctr(statei2.g, 1, 'r');
    ht = plot3Mat(plot_tip_traj, 'b', 1);
    ht2 = plot3Mat(plot_tip_traj2, 'r', 1);

    writeVideo(vid1, getframe(gcf));
    
end

close(vid1)


%% 2D plots

tip_traj = [];
tip_traj2 = [];
% tip_traj_est = [];
error_list = [];
error_list2 = [];
q_plot = [];
q_plot2 = [];

n_traj = length(plot_traj);

for ii = 1:n_traj
    stateii = plot_traj(ii);
    stateii2 = plot_traj2(ii);

    tip_traj = [tip_traj, stateii.p{3}(end, :)'];
    tip_traj2 = [tip_traj2, stateii2.p{3}(end, :)'];

    q_plot = [q_plot, stateii.q];
    q_plot2 = [q_plot2, stateii2.q];

    error_list = [error_list, stateii.error];
    error_list2 = [error_list2, stateii2.error];

end

% time serial
t_total = 10;

t0 = linspace(0, t_total, n_traj);

figure()
hold on
% plot(tip_traj_est(1, :))
% plot(tip_traj_est(2, :))
% plot(t0, 1000*x_path(2, :), 'g--', 'LineWidth',2);
% plot(t0, 1000*tip_traj(2, :), 'r')
plot(t0, 1000*x_path(2, idx_sampled), 'g--', 'LineWidth',2);
plot(t0, 1000*tip_traj(2, :), 'b')
plot(t0, 1000*tip_traj2(2, :), 'r')
xlabel('t (s)')
ylabel('Y (mm)')
% xlim([-0.2 4])
ylim([110 190])

plotConfig2D([600 300])
makeVid2D(gca, 'simu_out/vid_ctr_pid2_Y', [1 2])


figure()
hold on
plot(t0, q_plot([1 3 5], :)', 'b');
plot(t0, q_plot2([1 3 5], :)', 'r');

xlabel('t (s)')
ylabel('Y (mm)')

figure()
hold on
plot(s_list', 'r')
plot([state_traj.error], 'b')

% figure()
% hold on
% plot(t0, 1000*error_list)
% xlabel('t (s)')
% ylabel('Error (mm)')

% plotConfig2D([400 250])
% makeVid2D(gca, 'simu_out/vid_ctr_pid2_Y')

fprintf("Avg. tip error is %.3f mm, Max tip error is %.3f mm. \n", mean(1000*error_list), max(1000*error_list))

%% legend

figure()
hold on
legend box off
plot(nan, nan, 'b', 'LineWidth', 1)
plot(nan, nan, 'r', 'LineWidth', 1)
plot(nan, nan, 'g--', 'LineWidth',2)
legend('Cosserat-based control', 'LCP-based control', 'Desired trajectory')



%%
function q_sol = stepLCPControl(ctr, obj, x_desire, state_prev, d)

qc0 = state_prev.q;

    function y = obj_fun(q_135)

    q_full = [q_135(1), 0, q_135(2), 0, q_135(3), 0]';

    [state0, ~] = getFrictionalTipContactCTR(ctr, obj, q_full, state_prev, d);

    % y = norm(x_desire - state0.p{3}(end,:)') + 0.01*norm(q_135-qc0([1 3 5])' );
    y = norm(x_desire - state0.p{3}(end,:)');

    end

% q135_0 = [-4.7, -2.3, 0]';
% q135_0 = qc0([1 3 5]) + 0.01*[1 1 1]';

q135_0 = qc0([1 3 5]);

opt = optimset('Display','none', 'TolX', 1e-5, "TolFun", 1e-5);

AA = []; bb = []; 
AAeq = []; bbeq = []; 
% lb = [];
% ub = []; 
lb = q135_0 - 0.05*[1 1 1]';
ub = q135_0 + 0.05*[1 1 1]'; 

q135_sol = patternsearch(@obj_fun, q135_0, AA, bb, AAeq, bbeq, lb, ub, [], opt);
% q135_sol = fminsearch(@obj_fun, q135_0, opt);

% q135_0 = [-4.7, -2.3, 0]';
% 
% opt = optimoptions("fminunc", "Algorithm", "quasi-newton", "Display","iter", ...
%     "ObjectiveLimit", 1e-5, 'UseParallel',false);
% 
% q135_sol = fminunc(@obj_fun, q135_0, opt);

q_sol = [q135_sol(1), 0, q135_sol(2), 0, q135_sol(3), 0]';

end


function q_sol = stepLCPControl_mpc(ctr, obj, x_desire, state_prev0, d)

nh = size(x_desire, 2);   % horizon

    function y = obj_fun(q135_list_in)

        % qc = q135_list(:, 1);
        q135_list = reshape(q135_list_in, 3, nh);
        tip_traj = [];

        state_prev = state_prev0;

        for it = 1:nh

            qc135 = q135_list(:, it);
            qc_full = [qc135(1), 0, qc135(2), 0, qc135(3), 0]';

            [statec, ~] = getFrictionalTipContactCTR(ctr, obj, qc_full, state_prev, d);

            % update states
            state_prev = statec;

            tip_traj = [tip_traj, statec.p{3}(end, :)'];
        end

        y = norm(x_desire - tip_traj);

    end


q_ini = repmat(state_prev0.q([1 3 5]), [nh, 1]);
% q_ini = repmat([-4.7, -2.3, 0]', [1, nh]);

opt = optimset('Display','iter');

q135_sol = fminsearch(@obj_fun, q_ini, opt);

% opt = optimoptions("fminunc", "Algorithm", "quasi-newton", "Display","none", ...
%     "ObjectiveLimit", 1e-5);
% 
% q135_sol = fminunc(@obj_fun, q_ini, opt);

% excute first step
q_sol = [q135_sol(1,1), 0, q135_sol(2,1), 0, q135_sol(3,1), 0]';


end