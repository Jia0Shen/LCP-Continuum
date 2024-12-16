%% Test contact kinematics
clc;
clear all;
close all;

%% load python path for LCP solver
% pysys = py.sys.path;
% numpy_path = 'C:\Users\jshen359\AppData\Roaming\Python\Python38\site-packages';
% lcp_path = [pwd,'\LCP\lemkelcp'];
% pysys.append(numpy_path)
% pysys.append(lcp_path)
% np = py.importlib.import_module('numpy');
% lcp = py.importlib.import_module('lemkelcp');
%% Create tube
% n = 50;    % # of dicretization
% 
% L = 180;
% uhat = [0;1/70;0] * ones(1,n);
% uhat(:,1:n/2) = 0;

% tube = struct;
% tube.s = linspace(0, L, n);
% tube.uhat = uhat;
% 
% tube.kb = 1;
% tube.kt = 1/1.3;
% 
% tube.T_base = eye(4);

tube = CreatTube(150);

%% Create FRICTIONAL Obstacles
obstacles = cell(1,0);

capsules = cell(1,2);

p0 = [0;0;50];
r = 13;
height = 100;
z = [0;0;1];
mu = 0.5;
in_or_out = 'in';

R = LargeSO3([0,pi/2,0]);
p1 = p0 + (eye(3) + R) * [0;0;0.5*height];

capsules{1} = createCapsule(p0, z, r, height, mu, in_or_out);
capsules{2} = createCapsule(p1, R*z, r, height, mu, in_or_out);

% obstacles{1} = createAddedObject(capsules);
obstacles{1} = createPipe(capsules);

%% Create FRICTIONLESS Obstacles
obstacles_fl = obstacles;
for io = 1:length(obstacles_fl)
    obstacles_fl{io}.mu = 0.01;
end

%% Solve for shape without friction

[u, contacts, ~, R, p] = getInitialShape3(tube, obstacles);
% [T,R,p] = solveShape(tube.T_base, u, tube.s);

u_fl = u;
contacts_fl = contacts;
R_fl = R;
p_fl = p;



%% Plot
% obstacles2plot = obstacles;
obstacles2plot = capsules;

%test for contact detect
contactDetect = obstacles{1}.detectContact(tube, p, true);

figure('units','pixels','position',[0 0 840 840])
h = plot3Mat(p);
% set(h, 'LineWidth', 3);
hold on;
h_fl = plot3Mat(p_fl);
for io = 1:length(obstacles2plot)
    [X,Y,Z] = obstacles2plot{io}.getMesh(200);
    surf(X,Y,Z, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
    % surf(X,Y,Z, 'FaceAlpha', 0.5);
end
% h_dir = plot3Mat([zeros(3,1), [1;0;0]] * 10);

zlim([0 150])

axis equal
grid on

xlabel('x')
ylabel('y')

set(h, 'Color', 'r')
set(h_fl, 'Color', 'b')
% set(h, 'XData', p(1,:), 'YData', p(2,:), 'ZData', p(3,:));
% % set(h_fl, 'XData', p_fl(1,:), 'YData', p_fl(2,:), 'ZData', p_fl(3,:));
% dir = [tube.T_base(:,4), tube.T_base(:,4) + tube.T_base(:,1)*10];
% set(h_dir, 'XData', dir (1,:), 'YData', dir (2,:), 'ZData', dir (3,:), 'Marker','.');

view([0,-1,0])
axis([-60, 150, -120, 120, 0, 150])

drawnow



%% Solve for shape with friction

% creat video
vid = VideoWriter('fig_various_contact_vid.avi','Uncompressed AVI');
vid.FrameRate = 5;
open(vid)

tube.T_base(3,4) = 0;
num_int = 300;
beta_traj = [linspace(0,30,num_int), linspace(30,0,num_int)];
% remove redundant

beta_traj(diff(beta_traj) == 0) = [];
alpha_traj = zeros(size(beta_traj));

dalpha_traj = diff(alpha_traj);
dbeta_traj = diff(beta_traj);

t_total = 20;  % 8s

% the traj to record
tip_traj = p(:,end);
pc_traj = contacts.point;
nt = length(alpha_traj);

for i = 1:nt
    i;
    tic
    cornerFlag = true;
    MaxFricFlag = false;

    % record for debug.
    contact_prev = contacts;
    contact_prev_fl = contacts_fl;
    % u_prev = u;
    % p_prev = p;

    % ------ 1.Resolve contacts before we move ------ 
    contacts = [];
    [contacts, T, R, p] = detectAdditionalContacts(tube, u, obstacles, contacts, cornerFlag);

    if isempty(contacts)
        contacts = contact_prev;
    end

    contacts_fl = [];
    [contacts_fl, ~, R_fl, p_fl] = detectAdditionalContacts(tube, u_fl, obstacles, contacts_fl, cornerFlag);

    if isempty(contacts_fl)
        contacts_fl = contact_prev_fl;
    end

    % ------ 2. Move ------
    % tube.T_base(3,4) = i/10;

    alphaI = alpha_traj(i);
    betaI = beta_traj(i);

    tube.T_base = [RotZ(alphaI), [0,0,betaI]'; zeros(1,3), 1];
    if i < nt; tube.dalpha = dalpha_traj(i); tube.dbeta = dbeta_traj(i);
    else tube.dalpha = 0; tube.dbeta = 0; end

    % ------ 3. solve LCP and get contact force, p,u ------

    [u, contacts, ~, R, p] = getFrictionalContactShape3(tube, obstacles, u, contacts, cornerFlag);
    [u_fl, contacts_fl, ~, R_fl, p_fl] = getContactShape3(tube, obstacles, u_fl, contacts_fl, cornerFlag);
    % [u_fl, contacts_fl, ~, R_fl, p_fl] = getFrictionalContactShape3(tube,  obstacles_fl, u_fl, contacts_fl);

    % ------ 4. Output the solution at this step and plot ------
    if mod(i,1) == 0
        set(h, 'XData', p(1,:), 'YData', p(2,:), 'ZData', p(3,:));
        set(h_fl, 'XData', p_fl(1,:), 'YData', p_fl(2,:), 'ZData', p_fl(3,:));
        title('i='+string(i))

        % % draw the contact point and force
        % nc = length(contacts);
        % itc = zeros(1,nc);  % tube point ids on contact
        % p_contacts = zeros(3, nc);
        % f_contacts = zeros(3, nc);
        % n_contacts = zeros(3, nc);
        % for ic = 1:nc
        %     itc(ic) = contacts(ic).tube_point_id;
        %     p_contacts(:, ic) = contacts(ic).point;
        %     f_contacts(:, ic) = contacts(ic).force;
        %     n_contacts(:, ic) = contacts(ic).normal;
        % end

        drawnow

        writeVideo(vid,getframe(gcf));

        %warning if tip position has large error
        % if vecnorm(p_prev(:,end)-p(:,end)) > 10
        %     warning('Numerical Unstable at i=' +string(i))
        % end

    end

end

close(vid)

















