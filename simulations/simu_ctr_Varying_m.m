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
% obstacle{1} = createPlane(p0, n, mu);
% obstacle_fl{1} = createPlane(p0, n, 1e-5);

% tube.T_base(3,4) = 0;
q0 = [0 -0.1 -pi 0 0 0]';
% q1 = [0 0 pi 0 0 0]';
% q2 = [0 0 pi 0 2*pi 0]';
q1 = [0 -0.08 -pi 0 0 0]';
% q2 = [2*pi -0.08 pi 0 2*pi 0]';
q2 = q1 + (2*pi+pi/10)*[1 0 1 0 1 0]';
dalpha = deg2rad(1);   % rotate in 1 degree per step
dbeta = 2e-5;   % move in 1mm per step

q_traj1 = LinearInterplNdim([q0,q1], dbeta);
q_traj2 = LinearInterplNdim([q1,q2], dalpha);

q_traj = [q_traj1, q_traj2];

w0 = [0 0 0 0 0 0]';
yu0_guess = zeros(5,1);

%% 
n1 = size(q_traj1, 2);
n2 = size(q_traj2, 2);

n = size(q_traj, 2);

% d_list = [6 12];
d_list = [6 12 24 36 48 60 72 96 150];

tipTrajCell = {};
stateTrajCell = {};

% initialize
[~, p_ini, yu0, b_res, g_ini] = three_tube_fk(ctr, q0, w0, yu0_guess);

state_ini.p = p_ini;
state_ini.w = w0;
state_ini.yu0 = yu0;
state_ini.g = g_ini;
state_ini.q = q0;

for d = d_list

    tipContact = false;
    % yu0 = yu0_guess;

    tip_traj = [];
    state_traj = [];


    state_prev = state_ini;

    for i = 1:n

        qi = q_traj(:, i);
    
        if i > n1
            plane.mu = 0.3;
        else 
            plane.mu = mu;
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

    tipTrajCell{end+1} = tip_traj;
    stateTrajCell{end+1} = state_traj;

end

%% get traj and calculate accuracy
tipTraj_m6 = tipTrajCell{1};
tipTraj_m12 = tipTrajCell{2};
tipTraj_m150 = tipTrajCell{9};
stateTraj_m150 = stateTrajCell{9};

acc_list = [];

for i = 1:length(tipTrajCell)
    tipTraj_mi = tipTrajCell{i};

    % compare error with m=150
    acc = max(vecnorm(tipTraj_mi - tipTraj_m150, 2, 1));

    acc_list = [acc_list, acc];
end

figure()
plotConfig2D([400, 250])
plot(d_list, 1000*acc_list, '-r.', 'LineWidth', 1)
% plot(d_list, acc_list, '-r.', 'LineWidth', 1)
ylabel('Max Tip Error (mm)')
xlabel('Number of Edges')
ylim([0,max(1000*acc_list)*1.05])

saveas(gcf, "simu_out/fig_CTRCircle_Error", 'png')
close(gcf)

%%
figure()
plotConfig3D([400 400])
hold on
axis equal

pathCut = 1:n;
% pathCut = round(n/3):n;

plot(tipTraj_m150(1,pathCut), tipTraj_m150(2,pathCut), 'Color','r', 'LineWidth',1)
plot(tipTraj_m12(1,pathCut), tipTraj_m12(2,pathCut), 'Color','b', 'LineWidth',1)
plot(tipTraj_m6(1,pathCut), tipTraj_m6(2,pathCut), 'k', 'LineWidth',1)

xlim([-0.09, 0.09])
ylim([-0.09, 0.09])

saveas(gcf, "simu_out/fig_CTRCircle_tip.png", 'png')

makeVid2D(gca, "simu_out/vid_CTRCircle_tip")


%% plot the traj

vid1 = VideoWriter('simu_out/vid_ctr', 'MPEG-4');
vid1.FrameRate = n/10;
open(vid1);

figure()
hold on
plotConfig3D([420 640])

% plot3Mat(tip_traj)

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')
% axis equal
xlim([-0.1, 0.1])
ylim([-0.1, 0.1])
zlim([-0.02, 0.52])

view(45,10)

h = plotPlane(plane, 0.2);
h.FaceColor = 0.2*[1 1 1];
h.FaceAlpha = 0.2;

hs = []; ht1 = [];ht2 = [];ht3 = [];

for i = 1:n

    statei = stateTraj_m150(i);

    delete([hs ht1 ht2 ht3])

    % plot_three_tube(gi{:});
    hs = simple_plot_ctr(statei.g);
    ht1 = plot3Mat(tipTraj_m6(:,1:i), 'k');
    ht2 = plot3Mat(tipTraj_m12(:,1:i), 'b');
    ht3 = plot3Mat(tipTraj_m150(:,1:i), 'r');
    % ht3 = [];
    % pause(0.001)
    writeVideo(vid1, getframe(gcf));
    
end

close(vid1)

saveas(gcf, "simu_out/fig_CTRCircle.png", 'png')



%% legends
figure()
plotConfig2D([600 200])
hold on
plot(nan, nan, 'Color','r', 'LineWidth',1)
plot(nan, nan, 'b', 'LineWidth',1)
plot(nan, nan, 'k', 'LineWidth',1)
grid off
legend box off

legend('m=150', 'm=12', 'm=6', 'NumColumns', 3)
saveas(gcf, "simu_out/figs_CTRCircle_legend", 'png')
close(gcf)


