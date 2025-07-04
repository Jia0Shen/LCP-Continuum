clear; clc

ctr = build_ctr();

ctr_est = ctr;
ctr_est.E = ctr.E * (1 + 0.2);  % 5% uncertainty
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

q0 = [0 -0.1 -pi 0 0 0]';
q1 = [0 -0.085 -pi 0 0 0]';
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


%  def plane
wall_dis = p_ini{3}(end,3)+0.01;
p0 = [0;0;wall_dis];
n = [0;0;-1];
mu = 0.2;

plane = createPlane(p0, n, mu);

plane_fl = plane;
plane_fl.mu = 0;

% initialize
state_ini.p = p_ini;
state_ini.w = w0;
state_ini.yu0 = yu0;
state_ini.g = g_ini;
state_ini.q = q0;
state_prev = state_ini;

d = 60;

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

%% 

% design path
rp = norm(tip_traj_ini(1:2, end));
x_tip = tip_traj_ini(1, end);
y_tip = tip_traj_ini(2, end);
z_tip = tip_traj_ini(3, end);
theta0 = atan2(y_tip, x_tip);

n_path = 20;
directSign = 1;   % +1 for CCW, -1 for CW.
thetas = linspace(theta0, theta0+2*directSign*pi, n_path);
x_path = [rp*cos(thetas); rp*sin(thetas); z_tip*ones(1,n_path)];

% start the trajecotry
% save 
loadStates = false;
if loadStates
    state_load = load("simu_out\ctr_pdi_ini.mat");
    state_st = state_load.state_traj_ini(end);
else
    save("simu_out\ctr_pdi_ini.mat", "state_traj_ini");
    state_st = state_traj_ini(end);
end

%%

q_st = state_st.q;
w_st = state_st.w;
yu0_st = state_st.yu0;

xc = tip_traj_ini(:,end);
[J, C, ~, ~, B] = cal_Jacobian_ivp_one(ctr, q_st, w_st, yu0_st);
Jc = J(1:3,:);

% gain
K = 0*diag([1 1 1]);
nh = 100;   % horizon to move between waypoints.
d_points = norm(x_path(:,1) - x_path(:,2));
omega = 90;   % 90 deg/s
velo = 4*d_points*n_path/4;   % m/s
dt = 0.001;  % 1000 Hz
tipContact = true;
d = 60;
qc = q_st;
yu0_c = yu0_st;
state_c = state_st;
state_est = state_st;

state_traj = [];

% change the COF
% plane.mu = 0.15;
plane.mu = 0.1;

switchFlag = false;

pathIdx = [1:3,2,1,2];
for i = 1:length(pathIdx)   % [1:7, 6, 5, 4, 3, 4, 5, 4]

    % find jacobian - dq
    xd = x_path(:, pathIdx(i));

    if i == 4 || i == 6
        dt = 0;
        disp(length(state_traj));
    end

    for j = 1:nh

        % siwtch small only for 20 steps.
        if j > 10
            dt = 0.001;
        end

        % stop if reach the target.
        if norm(xd-xc) < 5e-4    % 0.5 mm accuracy
            fprintf('Reached the %d-th/%d target\n', i, length(pathIdx));
            break
        end

        %adaptive velocity. constant when far from the target.
        dist = norm(xd-xc);
        dist_tol = d_points/2;
        % dist_tol = d_points;
        if dist > dist_tol
            v0 = velo;
        else
            v0 = velo * dist/dist_tol;
            v0 = max(v0, velo/5);   % do not be too small
        end
        dx = v0*(xd-xc)/dist;   
        [J, C, ~, ~, B] = cal_Jacobian_ivp_one(ctr_est, state_est.q, state_est.w, state_est.yu0);
        % [J, C, ~, ~, B] = cal_Jacobian_ivp_one(ctr, state_c.q, w_st, yu0_c);
        
        % use all 6 DOF for j_idx = 1:6;
        j_idx = [1 3 5];
        % j_idx = 1:6;
        Jv = J(1:3,j_idx);

        % Jv_inv = Jv'*(inv(Jv*Jv' + 0.00001 * eye(3)));   % singularity-robust
        Jv_inv = Jv'*(inv(Jv*Jv' + 0.00001 * eye(3)));   % singularity-robust
        % Jv_inv = pinv(Jv);   % singularity-robust
        dq_idx = Jv_inv*(dx+K*(xd-xc));
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

        if tipContact
            % assume obs not moving here
            [state_c,tipContact] = getFrictionalTipContactCTR(ctr, plane, qc, state_prev, d);

            % estimate the state using frictionless model
            % [state_est,~] = getFrictionalTipContactCTR(ctr, plane, qc, state_est_prev, d);
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


end




%% plot the traj

% plot_traj = [state_traj_ini, state_traj];
% plot_traj = state_traj_ini;
plot_traj = state_traj;

plot_tip_traj = [];
plot_tip_traj_est = [];
n_traj = length(plot_traj);

% close(vid1)

vid1 = VideoWriter('simu_out/vid_ctr3', 'MPEG-4');
vid1.FrameRate = 110;
open(vid1);

figure()
hold on
plotConfig3D([640 640])

% plot3Mat(tip_traj)

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
axis equal
axislim = 0.08;
% xlim([-axislim, axislim])
% ylim([-axislim, axislim])
xlim([-axislim, 0.02])
ylim([-0.04, axislim])

zlim([-0.02, 0.52])


% view(30,45)
view(0, 90)

h = plotPlane(plane, 0.5);
h.FaceColor = 'k';
h.FaceAlpha = 0.2;
hs = []; ht = []; hs2 = []; ht2 = [];

% plot path
% plot3(x_path(1,:), x_path(2,:), x_path(3,:), 'bo')
plot3(x_path(1,1:7), x_path(2,1:7), x_path(3,1:7), 'ko')

for i = 1:n_traj

    statei = plot_traj(i);
    plot_tip_traj = [plot_tip_traj, statei.p{3}(end,:)'];
    plot_tip_traj_est = [plot_tip_traj_est, statei.g_est{3}(end, 13:15)'];

    delete([hs hs2 ht ht2])

    hs = simple_plot_ctr(statei.g, 1, 'r');
    % hs2 = simple_plot_ctr(statei.g_est);
    hs2 = simple_plot_ctr(statei.g_est, 1, 'b');
    % ht = plot3Mat(tip_traj_ini(:,1:i), 'b');
    ht = plot3Mat(plot_tip_traj, 'r');
    ht2 = plot3Mat(plot_tip_traj_est, 'b');

    writeVideo(vid1, getframe(gcf));
    
end

close(vid1)


%% 2D plots

tip_traj = [];
tip_traj_est = [];
q_plot = [];

for ii = 1:n_traj
    stateii = state_traj(ii);

    tip_traj = [tip_traj, stateii.p{3}(end, :)'];
    tip_traj_est = [tip_traj_est, stateii.g_est{3}(end, 13:15)'];
    q_plot = [q_plot, stateii.q];

end

figure()
hold on
% plot(tip_traj_est(1, :))
% plot(tip_traj_est(2, :))
plot(tip_traj(1, :), '.')
plot(tip_traj(2, :), '.')

figure()
hold on
plot(q_plot(1,:), '.')
plot(q_plot(3,:), '.')
plot(q_plot(5,:), '.')

%% legend

figure()
hold on
legend box off
plot(nan, nan, 'r', 'LineWidth', 1)
plot(nan, nan, 'b', 'LineWidth',1)
legend('Friction model', 'Frictionless model')
