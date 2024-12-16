%% Test contact kinematics
clc;
clear all;
% close all;

%% Create tube
% tube = CreatTube(150);

% outFBG = load('Sensing-FBGS\calib\reg_0417');
% curv0 = outFBG(end-209:end-1,1)/10;  % 1/cm -> 1/mm
% curvAngle0 = outFBG(end-209:end-1,2);
% s_vec = 0:208;
% R0 = eye(3); p0 = zeros(3,1);
% % reconstruct 3D precurvature
% [~, g0, uhat, g0_cos] = Curvature2shape(s_vec,curv0,curvAngle0, R0, p0);
% 
% % interplate to get more discretize points.
% s_m = linspace(0, 208, 400);
% uhat_m = interp1(s_vec, [uhat,zeros(3,1)]', s_m)';

% creat tube
% tube = CreatTube(200, s_vec, [uhat,zeros(3,1)]); % 213.88-14
% tube = CreatTube([], s_m, uhat_m);
% rotBase = 2.6;
tube = CreatTube(200);
rotBase = 0;
R0 = eye(3); p0 = zeros(3,1);

% physical parameters
tube.v = 0.3;  
tube.kb = 20.07e4;   % 20.07e-2 Nm2; 
tube.kt = tube.kb/(1+tube.v);   
% rotational compensation
PivotAngleX = - 0 /180*pi;
PivotAngleY = - 0 /180*pi;
DisplX = 0;
DisplY = 0;
tube.T_base = [RotX(PivotAngleX)*RotY(PivotAngleY)*RotZ(rotBase), p0+[DisplX;DisplY;0]; zeros(1,3), 1];

mu1 = 0.8; mu2 = 0.8;

CornerRange1 = 0.5;
CornerRange2 = 0.5;
MaxFricFlag = false;

%% Create FRICTIONAL Obstacles
obstacles = cell(1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = [0;0;1];
h1 = 1;
h2 = 110 - h1;
p1_cy = [0;0;170-h1/2];
p2_cy = [0;0;p1_cy(3)-(h1+h2)/2];
r = 4 - tube.rout;

% obstacles{end+1} = createCylinerCut(p0_cy, z, r, h, mu, 'in');
obstacles{1} = objCylinderCut(p1_cy, z, r, h1, mu1, 'in');
obstacles{2} = objCylinderCut(p2_cy, z, r, h2, mu2, 'in');
obstacles{2}.cornerFlag = false;

% tube.T_base(3,4) = 0;
num_int = 100;
angleLim = - pi/2;
insertLim = 0;
pn = 2;
alpha_traj = interp1(0:pn, -angleLim*rem(0:pn,2), linspace(0,pn,pn*num_int));
beta_traj = interp1(0:pn, -insertLim*rem(0:pn,2), linspace(0,pn,pn*num_int));
alpha_prev = alpha_traj(1);
t_total = 20;  % 8s

plotStep = 1;
pauseStep = 0.001;  % seconds
%% Solve for shape without friction

load_ini = false;
if ~load_ini
    [u, contacts, ~, R, p] = getInitialShapeIter(tube, obstacles, true);
    save('simu_out/simu_cylinder_ini.mat', "u", "contacts","R","p");  % ini3 for aligned axis
else
    cylinder_ini = load('simu_out/simu_cylinder_ini.mat', "u", "contacts","R","p");
    u = cylinder_ini.u;
    contacts = cylinder_ini.contacts;
    R = cylinder_ini.R;
    p = cylinder_ini.p;
end

% [u, contacts, ~, R, p] = getInitialShapeIter(tube, obstacles);

u_fl = u;
contacts_fl = contacts;
R_fl = R;
p_fl = p;

%% plot a sketch

obstacles2plot = obstacles;

figure()
hold on
% h = plot3Mat(p);
p_rot30 = RotZ(pi/6)*p;
p_rot30Minus = RotZ(-pi/6)*p;
h_rot1 = plot3Mat(p_rot30);
h_rot2 = plot3Mat(p_rot30Minus);

plotConfigCylinder([480 840], obstacles)

% set(h, 'Color', 'r')
set(h_rot1, 'Color', 'r', 'LineWidth',2)
set(h_rot2, 'Color', 'k', 'LineWidth',2)

view([-95, 55])
% camup([1 0 0])

axis([-10, 40, -30, 30, 30, 190])

saveas(gcf, "simu_out/fig_pipeSketch", 'png')

%% Plot
obstacles2plot = obstacles;

figure('units','pixels','position',[0 0 840 840])
hold on
% h_prev = plot3(nan,nan,nan, 'b');
h = plot3Mat(p);
h_fl = plot3Mat(p_fl);

plotConfigCylinder([480 840], obstacles)

set(h, 'Color', 'r')
set(h_fl, 'Color', 'k')
h_ars = plot3(nan,nan,nan, 'k');
hn_ars = plot3(nan,nan,nan, 'g');
h_c = plot3(nan,nan,nan,'k.');

view([0, 90])
axis([-10, 40, -10, 10, 0, 200])

drawnow

%% Solve for shape with friction

% creat video
vid = VideoWriter('simu_out/simu_cylinder_vid','MPEG-4');
vid.FrameRate = 5;
open(vid)

obstacles{1}.T_history(:,:,1) = [RotZ(alpha_traj(1)), p1_cy; 0 0 0 1];
obstacles{2}.T_history(:,:,1) = [RotZ(alpha_traj(1)), p2_cy; 0 0 0 1];

% the traj to record
tip_traj = p(:,end);
% contact_num = 2;

% if length(contacts) == 2
%     pc_traj = [contacts(1).point; contacts(2).point];
% end

pc_traj = [];
for ic = 1:length(contacts)
    pc_traj = [pc_traj; contacts(ic).point];
end

tip_traj_fl = p_fl(:,end);
% pc_traj_fl = contacts_fl.point;

nt = length(alpha_traj);

for i = 1:nt
    i;
    tic

    % record for debug.
    contact_prev = contacts;
    contact_prev_fl = contacts_fl;
    % u_prev = u;
    % p_prev = p;

    % ------ 1.Resolve contacts before we move ------ 

    % if i == 370
    %     disp('stop')
    % end

    contacts = [];
    [contacts, ~, ~, ~] = detectAdditionalContacts(tube, u, obstacles, contacts, CornerRange1);

    % contacts(1).force = contact_prev(1).force;

    if length(contacts) == 2
        disp("contact 1: " + contacts(1).type + ", pene: " + string(contacts(1).penetrateDepth) ...
            + ". contact 2: " + contacts(2).type + ", pene: " + string(contacts(2).penetrateDepth));
    end

    if ~strcmp(contacts(1).type, 'cornerContact')
        warning('contact type changed!')
    end

    % ------ 4. Output the solution at this step and plot ------
    if mod(i,plotStep) == 0
        % set(h_prev, 'XData', p(1,:), 'YData', p(2,:), 'ZData', p(3,:));
        % set(h_fl, 'XData', p_fl(1,:), 'YData', p_fl(2,:), 'ZData', p_fl(3,:));
        title('i='+string(i))
        % rotate the environment
        % set(h_env1, 'ZData', Z+betaI)
        % set(h_env2, 'ZData', Zline+betaI)
        % rotate([h_env1, h_env2], obstacles{1}.z, 180*(alphaI-alpha_prev)/pi)
        % alpha_prev = alphaI;

        % % draw the contact point and force
        nc = length(contacts);
        % itc = zeros(1,nc);  % tube point ids on contact
        p_contacts = zeros(3, nc);
        f_contacts = zeros(3, nc);
        n_contacts = zeros(3, nc);
        friction_ratio = zeros(1,nc);
        for ic = 1:nc
            % itc(ic) = contacts(ic).tube_point_id;
            fci = contact_prev(ic).force;
            nci = contact_prev(ic).normal;
            p_contacts(:, ic) = contacts(ic).point;
            f_contacts(:, ic) = fci;
            n_contacts(:, ic) = contacts(ic).normal;

            friction_ratio(ic) = norm((eye(3)-nci*nci')*fci) / norm(nci'*fci);
        end

        drawnow
        pause(pauseStep)

        writeVideo(vid,getframe(gcf));

        % draw arrow
        ar_start_c = p_contacts;
        ar_end_c = p_contacts + 8 * normalize(f_contacts, 'norm');
        ar_end_n = p_contacts + 4 * normalize(n_contacts, 'norm');

        set(h_c,'XData', p_contacts(1,:), 'YData', p_contacts(2,:), 'ZData', p_contacts(3,:));
        % set(h_ar,'XData', [ar_start_c(1,:)',ar_end_c(1,:)'], ...
        %          'YData', [ar_start_c(2,:)',ar_end_c(2,:)'], ...
        %          'ZData', [ar_start_c(3,:)',ar_end_c(3,:)']);
        for jj = 1:length(h_ars)
            delete(h_ars(jj));
            % delete(hn_ars(jj));
        end
        h_ars = plot3([ar_start_c(1,:);ar_end_c(1,:)],[ar_start_c(2,:);ar_end_c(2,:)],[ar_start_c(3,:);ar_end_c(3,:)],'g');
        % hn_ars = plot3([ar_start_c(1,:);ar_end_n(1,:)],[ar_start_c(2,:);ar_end_n(2,:)],[ar_start_c(3,:);ar_end_n(3,:)],'g');

        tip_traj = [tip_traj, p(:,end)];
        % tip_traj_fl = [tip_traj_fl, p_fl(:,end)];

        % if ~isempty(contacts)
        %     pc_traj = [pc_traj, contacts(1).point];
        %     % contacts.tube_point_id
        % end
        % 
        % if ~isempty(contacts_fl)
        %     pc_traj_fl = [pc_traj_fl, contacts_fl(1).point];
        %     % contacts_fl.tube_point_id
        % end

        % plot3(tip_traj(1,:), tip_traj(2,:), tip_traj(3,:), 'r')

        pc_traj = [pc_traj, reshape(p_contacts,[],1)];

        for ic = 1:length(contacts)
            % plot3(pc_traj(3*ic-2,:), pc_traj(3*ic-1,:),pc_traj(3*ic,:), 'c')
            % display the friction ratio.
            disp('Friction Ratio of contact ' + string(ic) + ' is ' + string(friction_ratio(ic)))
        end

     end

    % ------ 2. Move ------

    alphaI = alpha_traj(i);
    betaI = beta_traj(i);

    % tube.T_base = [RotZ(alphaI), [0,0,betaI]'; zeros(1,3), 1];

    % obstacles{1}.T_history(:,:,end+1) = [RotZ(alphaI), p0_cy; 0 0 0 1];
    obstacles{1}.T_history(:,:,end+1) = [RotZ(alphaI), p1_cy+[0;0;betaI]; 0 0 0 1];
    obstacles{1} = obstacles{1}.rebuild();

    obstacles{2}.T_history(:,:,end+1) = [RotZ(alphaI), p2_cy+[0;0;betaI]; 0 0 0 1];
    obstacles{2} = obstacles{2}.rebuild();

    % ------ 3. solve LCP and get contact force, p,u ------
    % obstacles{1}.mu = mu;
    [u, contacts, ~, R, p] = getFrictionalContactShape3(tube, obstacles, u, contacts, CornerRange2);
    % set mu = 0 for frictionless.
    % obstacles{1}.mu = 1e-2;
    % [u_fl, contacts_fl, ~, R_fl, p_fl] = getFrictionalContactShape3(tube, obstacles, u_fl, contacts_fl, CornerRange2);
    
    % [u_fl, contacts_fl, ~, R_fl, p_fl] = getContactShape3(tube,obstacles, u_fl, contacts_fl, CornerRange2);

    % add plot
    if mod(i,plotStep) == 0
        set(h, 'XData', p(1,:), 'YData', p(2,:), 'ZData', p(3,:));
        drawnow
    end

end

close(vid)

% ---------- draw the trajectory ----------
dt = t_total / length(alpha_traj);
% vc_traj = diff(pc_traj) / dt;
% vc_norm_traj = vecnorm(vc_traj, 2, 1);
vtip_traj = diff(tip_traj) / dt;
vtip_norm_traj = vecnorm(vtip_traj, 2, 1);

% 2.1 tip trajectory
figure()
hold on
plot(tip_traj(1,:), tip_traj(2,:), 'r')
% plot(tip_traj_fl(1,:), tip_traj_fl(2,:), 'b')
xlabel('x(mm)')
ylabel('y(mm)')
xlim([-0, 20])
ylim([-10, 10])
axis equal
title('tip trajectory')

figure();
hold on
plot(pc_traj(3,:))
plot(pc_traj(6,:))
plot(tip_traj(1,:))
plot(tip_traj(2,:))
ylabel('z(mm)')
xlabel('iter')

% % 2.2 contact point trajectory
% figure()
% hold on
% plot(pc_traj(1,:), pc_traj(2,:), 'r')
% plot(pc_traj_fl(1,:), pc_traj_fl(2,:), 'b')
% xlabel('x(mm)')
% ylabel('y(mm)')
% xlim([-120, 120])
% ylim([-120, 120])
% axis equal
% title('contact point trajectory')

% % 2.3. contact point x - time plot
% figure()
% hold on
% plot(dt*(0:nt), pc_traj(1,:), 'r')
% plot(dt*(0:nt), pc_traj_fl(1,:), 'b')
% xlabel('t(s)')
% ylabel('x(mm)')
% % xlim([-120, 120])
% ylim([-120, 120])
% title('contact point x-trajectory')
% 
% % 2.4. contact point y - time plot
% figure()
% hold on
% plot(dt*(0:nt), pc_traj(2,:), 'r')
% plot(dt*(0:nt), pc_traj_fl(2,:), 'b')
% xlabel('t(s)')
% ylabel('x(mm)')
% ylim([-120, 120])
% title('contact point y-trajectory')

% % 2.5. contact point x - alpha plot
% figure()
% hold on
% plot([0, alpha_traj], pc_traj(1,:), 'r')
% plot([0, alpha_traj], pc_traj_fl(1,:), 'b')
% xlabel('Base Rotation(rad)')
% ylabel('x(mm)')
% % xlim([-120, 120])
% ylim([-120, 120])
% title('contact point x-trajectory')
% 
% % 2.6. contact point y - alpha plot
% figure()
% hold on
% plot([0, alpha_traj], pc_traj(2,:), 'r')
% plot([0, alpha_traj], pc_traj_fl(2,:), 'b')
% xlabel('Base Rotation(rad)')
% ylabel('x(mm)')
% ylim([-120, 120])
% title('contact point y-trajectory')
% 
% figure()
% hold on
% plot([0, alpha_traj], atan2(pc_traj(2,:),pc_traj(1,:)), 'r')
% plot([0, alpha_traj], atan2(pc_traj_fl(2,:),pc_traj_fl(1,:)), 'b')
% xlabel('Base Rotation(rad)')
% ylabel('theta(rad)')
% ylim((angleLim+0.05)*[-1, 1])
% title('contact point theta-trajectory')

%% functions
function [u, contacts, T, R, p] = getInitialShapeIter(tube, obstacles, plotFlag)

if nargin == 2
    plotFlag = false;
end

% solve the shape from non-contact shape
h1 = obstacles{1}.h;  h2 = obstacles{2}.h;
p1_ini = [0; 0; 120-h1/2];
p2_ini = [0; 0; 120-h1-h2/2];
p1 = obstacles{1}.p0;
p2 = obstacles{2}.p0;
p1p2 = {p1_ini, p2_ini};

if p2(3) <= p2_ini(3)
    [u, contacts, T, R, p] = getInitialShapeEnergy(tube, obstacles);
else
    obstacles{1}.p0 = p1_ini;
    obstacles{2}.p0 = p2_ini;
    obstacles{1}.T_history(:,:,1) = [eye(3), p1_ini; 0 0 0 1];
    obstacles{2}.T_history(:,:,1) = [eye(3), p2_ini; 0 0 0 1];
    [u_ini, contacts_ini, ~, ~, ~] = getInitialShapeEnergy(tube, obstacles);
    beta_traj = 0:0.05:(p2(3)-p2_ini(3));
    nt = length(beta_traj);
    CornerRange1 = 0.5;
    CornerRange2 = 1;

    u = u_ini;    

if plotFlag
    figure('units','pixels','position',[0 0 840 840])
    h = plot3(nan,nan,nan, 'r');
    hold on;
    [X1,Y1,Z1] = obstacles{1}.getMesh(50);
    [X2,Y2,Z2] = obstacles{2}.getMesh(50);
    h_env1 = surf(X1,Y1,Z1, 'FaceColor', 'k',  'FaceAlpha', 0.3, 'EdgeColor', 'k', 'EdgeAlpha',0.6, 'LineWidth',0.5, 'MeshStyle','row');
    h_env2 = surf(X2,Y2,Z2, 'FaceColor', 'k',  'FaceAlpha', 0.3, 'EdgeColor', 'k', 'EdgeAlpha',0.6, 'LineWidth',0.5, 'MeshStyle','row');

    
    axis equal
    grid on
    xlabel('x')
    ylabel('y')

    h_ars = plot3(nan,nan,nan, 'k');
    hn_ars = plot3(nan,nan,nan, 'g');
    h_c = plot3(nan,nan,nan,'k.');
    
    view([0, 10])
    axis([-25, 25, -25, 25, 0, 200])
    
    drawnow
end

    for i = 1:nt
    
        contacts = [];
        [contacts, ~, ~, ~] = detectAdditionalContacts(tube, u, obstacles, contacts, CornerRange1);

        % if length(contacts) == 1
        %     warning('2 contact')
        %     % disp(contacts.type)
        % end

        for io = 1:2
            pi_ini = p1p2{io};
            obstacles{io}.T_history(:,:,end+1) = [eye(3), pi_ini+[0;0;beta_traj(i)]; 0 0 0 1];
            obstacles{io} = obstacles{io}.rebuild();
            obstacles{io}.mu = 0.1;
        end
    
        % obstacles{1}.mu = 0.01;
        [u, contacts, T, R, p] = getFrictionalContactShape3(tube, obstacles, u, contacts, CornerRange2);

        if plotFlag
            set(h, 'XData', p(1,:), 'YData', p(2,:), 'ZData', p(3,:));
            title('i='+string(i))
            set(h_env1, 'ZData', Z1+beta_traj(i))
            set(h_env2, 'ZData', Z2+beta_traj(i))
            drawnow

            % draw contacts.
            nc = length(contacts);
            p_contacts = zeros(3, nc);
            f_contacts = zeros(3, nc);
            n_contacts = zeros(3, nc);
            for ic = 1:nc
                p_contacts(:, ic) = contacts(ic).point;
                f_contacts(:, ic) = contacts(ic).force;
                n_contacts(:, ic) = contacts(ic).normal;
            end

            % draw arrow
            ar_start_c = p_contacts;
            ar_end_c = p_contacts + 15 * normalize(f_contacts, 'norm');
            ar_end_n = p_contacts + 7.5 * normalize(n_contacts, 'norm');

            set(h_c,'XData', p_contacts(1,:), 'YData', p_contacts(2,:), 'ZData', p_contacts(3,:));
            % set(h_ar,'XData', [ar_start_c(1,:)',ar_end_c(1,:)'], ...
            %          'YData', [ar_start_c(2,:)',ar_end_c(2,:)'], ...
            %          'ZData', [ar_start_c(3,:)',ar_end_c(3,:)']);
            for jj = 1:length(h_ars)
                delete(h_ars(jj));
                delete(hn_ars(jj));
            end
            h_ars = plot3([ar_start_c(1,:);ar_end_c(1,:)],[ar_start_c(2,:);ar_end_c(2,:)],[ar_start_c(3,:);ar_end_c(3,:)],'r');
            hn_ars = plot3([ar_start_c(1,:);ar_end_n(1,:)],[ar_start_c(2,:);ar_end_n(2,:)],[ar_start_c(3,:);ar_end_n(3,:)],'g');


        end
    
    end

end

end