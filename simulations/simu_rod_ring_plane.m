%% Test contact kinematics
clc;
clear all;
close all;

%% Create tube
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

mu1 = 0.5; mu2 = 0.2;  % mu1 for corner; mu2 for plane
    
CornerRange1 = 0.5;
CornerRange2 = 0.5;
MaxFricFlag = false;

%% Create FRICTIONAL Obstacles
obstacles = cell(1,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z = [0;0;1];
h1 = 3;
p1_cy = [0;0;160.1-h1/2];
r = 10;

obstacles{1} = objCylinderCut(p1_cy, z, r, h1, mu1, 'in');
obstacles{1}.cornerFlag = true;

wall_distance = 170;
p0_plane = [0;0;wall_distance];
z_plane = [0;0;-1];

obstacles{2} = objPlane(p0_plane, z_plane, mu2);

% tube.T_base(3,4) = 0;
num_int = 100;
angleLim = pi/2;
insertLim = 0;
pn = 2;
alpha_traj = interp1(0:pn, -angleLim*rem(0:pn,2), linspace(0,pn,pn*num_int));
beta_traj = interp1(0:pn, -insertLim*rem(0:pn,2), linspace(0,pn,pn*num_int));
alpha_prev = alpha_traj(1);
t_total = 20;  % 8s

plotStep = 1;
pauseStep = 0.001;  % seconds
load_ini = true;
%% Solve for shape without friction

if ~load_ini
    [u, contacts, ~, R, p] = getInitialShapeIter(tube, obstacles, true);
    save('simu_out/simu_cornerPlane_ini.mat', "u", "contacts","R","p");  % ini3 for aligned axis
else
    cylinder_ini = load('simu_out/simu_cornerPlane_ini.mat', "u", "contacts","R","p");
    u = cylinder_ini.u;
    contacts = cylinder_ini.contacts;
    R = cylinder_ini.R;
    p = cylinder_ini.p;
end

u_fl = u;
contacts_fl = contacts;
R_fl = R;
p_fl = p;


%% Plot
% obstacles2plot = obstacles;
obstacles2plot = obstacles;

rp = 60;
figure('units','pixels','position',[0 0 840 840])
hold on
% h_prev = plot3(nan,nan,nan, 'b');
h = plot3Mat(p);
h_fl = plot3Mat(p_fl);

[X1,Y1,Z1] = obstacles2plot{1}.getMesh(50);
[X2,Y2,Z2] = obstacles2plot{2}.getMesh(rp,100);
h_env1 = surf(X1,Y1,Z1, 'FaceColor', 'k',  'FaceAlpha', 0.3, 'EdgeColor', 'k', 'MeshStyle','row', 'EdgeAlpha',0.6, 'LineWidth',0.5);
h_env2 = surf(X2,Y2,Z2, 'FaceColor', 'k',  'FaceAlpha', 0.3, 'EdgeColor', 'k', 'MeshStyle','row', 'EdgeAlpha',0.6, 'LineWidth',0.5);
fill3(X2(1,:),Y2(1,:),Z2(1,:),'k','FaceAlpha', 0.3)
fill3(X2(2,:),Y2(2,:),Z2(2,:),'k','FaceAlpha', 0.3)
h_pin = plot3([p0_plane(1),rp],[p0_plane(2),0],[p0_plane(3)+5,p0_plane(3)+5], 'k--', 'LineWidth', 2);

axis equal
grid on
xlabel('x')
ylabel('y')

set(h, 'Color', 'r')
set(h_fl, 'Color', 'k')
h_ars = plot3(nan,nan,nan, 'k');
hn_ars = plot3(nan,nan,nan, 'g');
h_c = plot3(nan,nan,nan,'k.');

view([10, 50])
axis([-70, 70, -70, 70, 0, 200])

drawnow

%% Solve for shape with friction

% creat video
vid = VideoWriter('simu_out/simu_cornerPlane_vid','MPEG-4');
vid.FrameRate = 5;
open(vid)

obstacles{1}.T_history(:,:,1) = [RotZ(alpha_traj(1)), p1_cy; 0 0 0 1];
obstacles{2}.T_history(:,:,1) = [RotZ(alpha_traj(1)), p0_plane; 0 0 0 1];

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
FricRatioPath = [0, 0];

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

    % if isempty(contacts)
    %     contacts = contact_prev;
    % end

    % if length(contacts) == 2
    %     disp("contact 1: " + contacts(1).type + ", pene: " + string(contacts(1).penetrateDepth) ...
    %         + ". contact 2: " + contacts(2).type + ", pene: " + string(contacts(2).penetrateDepth));
    % end

    % if ~strcmp(contacts(1).type, 'cornerContact')
    %     warning('contact type changed!')
    % end

    % ------ 4. Output the solution at this step and plot ------
    if mod(i,plotStep) == 0
        % set(h_prev, 'XData', p(1,:), 'YData', p(2,:), 'ZData', p(3,:));
        % set(h_fl, 'XData', p_fl(1,:), 'YData', p_fl(2,:), 'ZData', p_fl(3,:));
        title('i='+string(i))
        % rotate the environment
        % set(h_env1, 'ZData', Z+betaI)
        % set(h_env2, 'ZData', Zline+betaI)       

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

        FricRatioPath = [FricRatioPath; friction_ratio];

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
            delete(hn_ars(jj));
        end
        h_ars = plot3([ar_start_c(1,:);ar_end_c(1,:)],[ar_start_c(2,:);ar_end_c(2,:)],[ar_start_c(3,:);ar_end_c(3,:)],'c');
        hn_ars = plot3([ar_start_c(1,:);ar_end_n(1,:)],[ar_start_c(2,:);ar_end_n(2,:)],[ar_start_c(3,:);ar_end_n(3,:)],'g');

        tip_traj = [tip_traj, p(:,end)];
        % tip_traj_fl = [tip_traj_fl, p_fl(:,end)];

        plot3(tip_traj(1,:), tip_traj(2,:), tip_traj(3,:), 'r')

        pc_traj = [pc_traj, reshape(p_contacts,[],1)];

        for ic = 1:length(contacts)
            plot3(pc_traj(3*ic-2,:), pc_traj(3*ic-1,:),pc_traj(3*ic,:), 'c')
            % display the friction ratio.
            disp('Friction Ratio of contact ' + string(ic) + ' is ' + string(friction_ratio(ic)))
        end

     end

    % ------ 2. Move ------

    alphaI = alpha_traj(i);
    betaI = beta_traj(i);

    obstacles{2}.T_history(:,:,end+1) = [RotZ(alphaI), p0_plane+[0;0;betaI]; 0 0 0 1];
    obstacles{2} = obstacles{2}.rebuild();

    % ------ 3. solve LCP and get contact force, p,u ------
    % obstacles{1}.mu = mu;
    [u, contacts, ~, R, p] = getFrictionalContactShape3(tube, obstacles, u, contacts, CornerRange2);

    % add plot
    if mod(i,plotStep) == 0
        set(h, 'XData', p(1,:), 'YData', p(2,:), 'ZData', p(3,:));
        rotate(h_env2, [0;0;1], 180*(alphaI-alpha_prev)/pi, p1_cy')
        rotate(h_pin, [0;0;1], 180*(alphaI-alpha_prev)/pi, p1_cy')
        alpha_prev = alphaI;
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

% 2.2. friction ratio path
% figure()
% hold on
% plot(dt*(0:nt), FricRatioPath(:,1), 'r')
% plot(dt*(0:nt), FricRatioPath(:,2), 'b')
% ylim([0,1])
% xlabel('t(s)')
% ylabel('\mu')

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
p_plane = obstacles{2}.p0;
pplane_ini = [0;0;185];
p_cy = obstacles{1}.p0;
pcy_ini = [0;0;140];

u = tube.uhat;

movePlaneFlag = true;  moveCylinderFlag = true;
if p_plane(3) >= pplane_ini(3)
    movePlaneFlag = false;
    pplane_ini = pplane;
end

if p_cy(3) <= pcy_ini(3)
    moveCylinderFlag = false;
    pcy_ini = p_cy;
end

% plan the insertion motion
insertStepLen = 0.1;
insertStepNum = max(pplane_ini(3)-p_plane(3), p_cy(3)-pcy_ini(3)) / insertStepLen;
betaCy = linspace(0, p_cy(3)-pcy_ini(3), insertStepNum);
betaPlane = linspace(0, p_plane(3)-pplane_ini(3), insertStepNum);

if ~movePlaneFlag && ~moveCylinderFlag
    [contacts, T, R, p] = detectAdditionalContacts(tube, u, obstacles, [], 1);
    if ~isempty(contacts)
        error('inital shape detects contacts, change p0_ini!')
    end
else
    % use initial p0 to iteratively solve the shape
    obstacles{1}.p0 = pcy_ini;
    obstacles{2}.p0 = pplane_ini;
    [contacts, T, R, p] = detectAdditionalContacts(tube, u, obstacles, [], 1);
    if ~isempty(contacts)
        error('inital shape detects contacts, change p0_ini!')
    end
    obstacles{1}.T_history(:,:,1) = [eye(3), pcy_ini; 0 0 0 1];
    obstacles{2}.T_history(:,:,1) = [eye(3), pplane_ini; 0 0 0 1];
    CornerRange1 = 0.5;
    CornerRange2 = 0.5;

    if plotFlag
        figure('units','pixels','position',[0 0 840 840])
        h = plot3(nan,nan,nan, 'r');
        hold on;
        [X1,Y1,Z1] = obstacles{1}.getMesh(50);
        [X2,Y2,Z2] = obstacles{2}.getMesh(60, 50);
        h_env1 = surf(X1,Y1,Z1, 'FaceColor', 'k',  'FaceAlpha', 0.3, 'EdgeColor', 'k', 'EdgeAlpha',0.6, 'LineWidth',0.5, 'MeshStyle','row');
        h_env2 = surf(X2,Y2,Z2, 'FaceColor', 'k',  'FaceAlpha', 0.3, 'EdgeColor', 'k', 'EdgeAlpha',0.6, 'LineWidth',0.5, 'MeshStyle','row');
        patchButt = fill3(X2(1,:),Y2(1,:),Z2(1,:),'k','FaceAlpha', 0.3);
        patchTop = fill3(X2(2,:),Y2(2,:),Z2(2,:),'k','FaceAlpha', 0.3);

        axis equal
        grid on
        xlabel('x')
        ylabel('y')
    
        h_ars = plot3(nan,nan,nan, 'k');
        hn_ars = plot3(nan,nan,nan, 'g');
        h_c = plot3(nan,nan,nan,'k.');
        
        view([0, 10])
        axis([-25, 75, -25, 25, 0, 200])
        
        drawnow
    end

    for i = 1:insertStepNum
    
        contacts = [];
        [contacts, ~, ~, ~] = detectAdditionalContacts(tube, u, obstacles, contacts, CornerRange1);

        % if length(contacts) == 1
        %     warning('2 contact')
        %     % disp(contacts.type)
        % end

        % if i == 302
        %     disp('')
        % end

        obstacles{1}.T_history(:,:,end+1) = [eye(3), pcy_ini+[0;0;betaCy(i)]; 0 0 0 1];
        obstacles{1} = obstacles{1}.rebuild();
        obstacles{1}.mu = 0.01;

        obstacles{2}.T_history(:,:,end+1) = [eye(3), pplane_ini+[0;0;betaPlane(i)]; 0 0 0 1];
        obstacles{2} = obstacles{2}.rebuild();
        obstacles{2}.mu = 0.01;
        
        [u, contacts, T, R, p] = getFrictionalContactShape3(tube, obstacles, u, contacts, CornerRange2);

        if plotFlag
            set(h, 'XData', p(1,:), 'YData', p(2,:), 'ZData', p(3,:));
            title('i='+string(i))
            set(h_env1, 'ZData', Z1+betaCy(i))
            set(h_env2, 'ZData', Z2+betaPlane(i))
            set(patchButt, 'ZData', Z2(1,:)+betaPlane(i))
            set(patchTop, 'ZData', Z2(2,:)+betaPlane(i))
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