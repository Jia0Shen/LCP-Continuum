clear; clc

ctr = build_ctr();

%% Create Obstacles, all in meter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

wall_dis = 0.5;
p0 = [0;0;wall_dis];
n = [0;0;-1];
mu = 0.1;

plane = createPlane(p0, n, mu);
plane_fl = createPlane(p0, n, 1e-5);

q0 = [0 -0.1 -pi 0 0 0]';
q1 = [0 -0.08 -pi 0 0 0]';
q2 = [2*pi -0.08 pi 0 2*pi 0]';
dalpha = deg2rad(1);   % rotate in 1 degree per step
dbeta = 2e-5;   % move in 2e-2 mm per step

q_traj1 = LinearInterplNdim([q0,q1], dbeta);
q_traj2 = LinearInterplNdim([q1,q2], dalpha);

q_traj = [q_traj1, q_traj2];

w0 = [0 0 0 0 0 0]';
yu0_guess = zeros(5,1);

%% 
n1 = size(q_traj1, 2);
n2 = size(q_traj2, 2);

n = size(q_traj, 2);

%start 
tic

tipContact = false;
% yu0 = yu0_guess;

tip_traj = [];
tip_traj_fl = [];

state_traj = [];
state_traj_fl = [];

% initialize
[~, p_ini, yu0, b_res, g_ini] = three_tube_fk(ctr, q0, w0, yu0_guess);

state_ini.p = p_ini;
state_ini.w = w0;
state_ini.yu0 = yu0;
state_ini.g = g_ini;
state_ini.q = q0;
state_prev = state_ini;

d = 60;

for i = 1:n

    qi = q_traj(:, i);

    if i > n1
        plane.mu = 0.3;
    end

    if tipContact
        % assume obs not moving here
        [state_update,tipContact] = getFrictionalTipContactCTR(ctr, plane, qi, state_prev, d);

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

    tip_traj = [tip_traj, p_tip];
    state_traj = [state_traj, state_update];

end

tolTime = toc;

disp(n/tolTime)

%% plot the traj

vid1 = VideoWriter('simulations/simu_out/vid_ctr', 'MPEG-4');
open(vid1);

figure()
hold on
plotConfig3D([640 640])

% plot3Mat(tip_traj)

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
xlim([-0.15, 0.15])
ylim([-0.15, 0.15])
zlim([-0.02, 0.6])

view(45,10)

h = plotPlane(plane, 0.5);
h.FaceColor = 'g';
h.FaceAlpha = 0.2;
hs = []; ht = [];

for i = 1:n

    statei = state_traj(i);

    delete([hs ht])

    hs = simple_plot_ctr(statei.g);
    ht = plot3Mat(tip_traj(:,1:i), 'b');

    writeVideo(vid1, getframe(gcf));
    
end

close(vid1)


%% 2D plots


