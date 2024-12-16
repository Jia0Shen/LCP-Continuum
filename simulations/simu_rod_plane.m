clear; clc
close all

tube = CreatTube(150);

R_base = [0 0 1; 0 -1 0; 1 0 0];
tube.T_base = [R_base, zeros(3,1); zeros(1,3), 1];

obstacles = cell(1,0);

wall_dis = 10;
p0 = [0;0;wall_dis];
n = [0;0;-1];
mu = 2.8;

obstacles{1} = createPlane(p0, n, mu);

% solve for shape without friction.

[u, contacts, ~, R, p] = getInitialShape3(tube, obstacles);
% [T,R,p] = solveShape(tube.T_base, u, tube.s);

% solve the LCP to correct
CornerRange1 = 0.5;
CornerRange2 = 0.5;

contacts = [];
[contacts, ~, ~, ~] = detectAdditionalContacts(tube, u, obstacles, contacts, CornerRange1);

obstacles{1}.mu = 0.001;
[u, contacts, T, R, p] = getFrictionalContactShape3(tube, obstacles, u, contacts, CornerRange2);
obstacles{1}.mu = mu;

contacts_fl = contacts;
R_fl = R;
p_fl = p;


%% Plot
% obstacles2plot = obstacles{1};
% 
% figure('units','pixels','position',[0 0 600 600])
% h = plot3Mat(p);
% hold on;
% h_fl = plot3Mat(p_fl);
% for io = 1:length(obstacles2plot)
%     [X,Y,Z] = obstacles2plot.getMesh(300);
%     surf(X,Y,Z, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
% end
% % h_dir = plot3Mat([zeros(3,1), [1;0;0]] * 10);
% 
% axis equal
% grid on
% 
% xlabel('x')
% ylabel('y')
% 
% set(h, 'Color', 'r')
% set(h_fl, 'Color', 'b')
% % set(h, 'XData', p(1,:), 'YData', p(2,:), 'ZData', p(3,:));
% % % set(h_fl, 'XData', p_fl(1,:), 'YData', p_fl(2,:), 'ZData', p_fl(3,:));
% % dir = [tube.T_base(:,4), tube.T_base(:,4) + tube.T_base(:,1)*10];
% % set(h_dir, 'XData', dir (1,:), 'YData', dir (2,:), 'ZData', dir (3,:), 'Marker','.');
% 
% view([0,45])
% axis([-120, 120, -120, 120, 0, wall_dis+5])

%% move the base and solve for the shape

nt = 1000;
beta_traj = linspace(0,20, nt);

beta0 = beta_traj(1);

% ini
tip_traj = zeros(3,nt);
whole_traj = zeros(3*size(p,2), nt);
contact_traj = cell(1,nt);

tip_traj(:,1) = p(:,end);
whole_traj(:,1) = reshape(p,[],1);

for i = 2:nt

    contacts = [];
    [contacts, ~, ~, ~] = detectAdditionalContacts(tube, u, obstacles, contacts, CornerRange1);
    contact_traj{i-1} = contacts;

    % actuation along X positive
    betaI = beta_traj(i);
    tube.T_base(1,4) = betaI;    

    [u, contacts, T, R, p] = getFrictionalContactShape3(tube, obstacles, u, contacts, CornerRange2);

    tip_traj(:,i) = p(:,end);
    whole_traj(:,i) = reshape(p,[],1);

end

%% plot the traj

obstacles2plot = obstacles{1};

figure('units','pixels','position',[0 0 700 400])
plotConfig3D([700 400])
hold on;
axis equal
grid on
box on

xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

% set(h, 'XData', p(1,:), 'YData', p(2,:), 'ZData', p(3,:));
% % set(h_fl, 'XData', p_fl(1,:), 'YData', p_fl(2,:), 'ZData', p_fl(3,:));
% dir = [tube.T_base(:,4), tube.T_base(:,4) + tube.T_base(:,1)*10];
% set(h_dir, 'XData', dir (1,:), 'YData', dir (2,:), 'ZData', dir (3,:), 'Marker','.');
% plot3(dir(1,:),dir(2,:),dir(3,:))
view([0,0])
axis([-10, 160, -120, 120, -30, wall_dis+15])

% loop for the trajectory
color = [0 0.4470 0.7410];
firstFlag = true;

idx_trans = 346;   % HARD CODE!
p_trans = reshape(whole_traj(:,idx_trans), 3, []);

for i = [10:100:idx_trans-1,idx_trans,idx_trans+1:100:nt]
    p_i = reshape(whole_traj(:,i), 3, []);
    tip_x = tip_traj(1,i);

    if i >= idx_trans
        color = [1 0 0];
    end
    % if abs(p_fl(1,end)-tip_x) > 0.01 && firstFlag
    %     % transform to sliding
    %     color = [1 0 0];
    %     idx_trans = i;
    %     p_trans = p_i;
    %     firstFlag = false;
    % end

    plot3Mat(p_i, color)
end

p_end = reshape(whole_traj(:,end), 3, []);

plot3Mat(p_fl, 'k', 2);
plot3Mat(p_trans, [0 0.4470 0.7410], 2);
plot3Mat(p_end, 'r', 2);

% draw the base
ar = 5;
dir_fl = [p_fl(:,1), p_fl(:,1) - ar*n];
dir_trans = [p_trans(:,1), p_trans(:,1) - ar*n];
dir_end = [p_end(:,1), p_end(:,1) - ar*n];

plot3Mat(dir_fl, 'k', 0.25);
plot3Mat(dir_trans, [0 0.4470 0.7410], 0.25);
plot3Mat(dir_end, 'r', 0.25);

% saveas(gcf, "simu_outs/figs_case1/fig_case1", 'png')

%% make video 

figure('units','pixels','position',[0 0 700 400])
plotConfig3D([700 400])
hold on;
axis equal
grid on
box on

xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

view([0,0])
axis([-10, 160, -120, 120, -30, wall_dis+15])

totalTime = 10;  % 10s
vidpath = "simu_out/vid_rod_plane";
vid = VideoWriter(vidpath,'MPEG-4');
vid.FrameRate = round(nt/totalTime);
open(vid)

pi = plot3(nan,nan,nan,'r', 'LineWidth',2);
Dir = plot3(nan,nan,nan,'r', 'LineWidth',1);

for i = 1:nt
    p_i = reshape(whole_traj(:,i), 3, []);

    % arrow at the base
    ar = 5;
    dir_trans = [p_i(:,1), p_i(:,1) - ar*n];

    set(pi, 'XData', p_i(1,:), ...
        'YData', p_i(2,:), ...
        'ZData', p_i(3,:));

    set(Dir, 'XData', dir_trans(1,:), ...
        'YData', dir_trans(2,:), ...
        'ZData', dir_trans(3,:));

    drawnow
    writeVideo(vid,getframe(gcf));
end

% draw the base
% ar = 5;
% dir_fl = [p_fl(:,1), p_fl(:,1) - ar*n];
% dir_trans = [p_trans(:,1), p_trans(:,1) - ar*n];
% dir_end = [p_end(:,1), p_end(:,1) - ar*n];
% 
% plot3Mat(dir_fl, 'k', 0.25);
% plot3Mat(dir_trans, [0 0.4470 0.7410], 0.25);
% plot3Mat(dir_end, 'r', 0.25);

close(vid)
%% 
idx_trans = 346;   % HARD CODE!

figure()
plotConfig2D([700, 200]);
plot(beta_traj(1:idx_trans), tip_traj(1,1:idx_trans), 'Color', [0 0.4470 0.7410], 'LineWidth', 1)
plot(beta_traj(idx_trans:end), tip_traj(1,idx_trans:end), 'r', 'LineWidth', 1)
plot(beta_traj(idx_trans), tip_traj(1,idx_trans), 'r.')
xlabel('Insertion (mm)')
ylabel('x (mm)')
ylim([135, 155])
xlim([0, 21])

% saveas(gcf, "simu_outs/figs_case1/fig_case1_traj", 'png')

%% make video
figure()
plotConfig2D([700, 200]);
plot(beta_traj, tip_traj(1,:), 'Color', [1 0 0], 'LineWidth', 1)
xlabel('Insertion (mm)')
ylabel('x (mm)')
ylim([135, 155])
xlim([0, 21])

% makeVid2D(gca, "simu_outs/figs_case1/vid_case1_traj")

%% legend
% figure()
% plotConfig2D([200, 200]);
% box off
% grid off
% plot(nan,nan,'Color', 'k', 'LineWidth',2);
% plot(nan,nan,'Color', [0 0.4470 0.7410], 'LineWidth',2);
% plot(nan,nan,'Color', 'r', 'LineWidth',2);
% 
% legend('Initial shape', 'Static friction',  'Sliding friction')
% 
% % saveas(gcf, "simu_outs/figs_case1/fig_case1_legend", 'png')