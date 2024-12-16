clear; clc
close all

tube = CreatTube(150);

R_base = [0 0 1; 0 -1 0; 1 0 0];
tube.T_base = [R_base, zeros(3,1); zeros(1,3), 1];

obstacles = cell(1,0);
obstacles_fl = cell(1,0);

wall_dis = 10;
p0 = [0;0;wall_dis];
n = [0;0;-1];
mu = 0.5;

obstacles{1} = createPlane(p0, n, mu);
obstacles_fl{1} = createPlane(p0, n, 1e-5);

% solve for shape without friction.

[u, contacts, ~, R, p] = getInitialShape3(tube, obstacles_fl);
% [T,R,p] = solveShape(tube.T_base, u, tube.s);

% solve the LCP to correct
CornerRange1 = 0.5;
CornerRange2 = 0.5;

contacts = [];
[contacts, ~, ~, ~] = detectAdditionalContacts(tube, u, obstacles_fl, contacts, CornerRange1);

[u, contacts, T, R, p] = getFrictionalContactShape3(tube, obstacles_fl, u, contacts, CornerRange2);


contacts_fl = contacts;
R_fl = R;
p_fl = p;
u_fl = u;

%%

nt = 200;
alpha_traj = linspace(0,3*pi, nt);

R_cir = 50;
x_traj = R_cir * sin(alpha_traj);
y_traj = R_cir * cos(alpha_traj);
act_traj = [x_traj; y_traj; zeros(size(x_traj))];

alpha0 = alpha_traj(1);

% param.d = 150;
d_list = 6:6:150;

% try different number of friction cone edges
results_mat = cell(1,length(d_list));

for it = 1:length(d_list)
    param.d = d_list(it);

    % param.d = 150;
    tic

    % ini
    tip_trajI = zeros(3,nt);
    whole_trajI = zeros(3*size(p,2), nt);
    contact_trajI = cell(1,nt);

    tip_trajI(:,1) = p(:,end);
    whole_trajI(:,1) = reshape(p,[],1);

    for i = 2:nt

        contacts = [];
        [contacts, ~, ~, ~] = detectAdditionalContacts(tube, u, obstacles, contacts, CornerRange1);

        % actuation along circular path
        alphaI = alpha_traj(i);
        tube.T_base(1:2,4) = [x_traj(i), y_traj(i)]';  
    
        [u, contacts, ~, ~, p] = getFrictionalContactShape3(tube, obstacles, u, contacts, CornerRange2, param);

        contact_trajI{i-1} = contacts;

        tip_trajI(:,i) = p(:,end);
        whole_trajI(:,i) = reshape(p,[],1);
    
    end


    TotalSpeedI = toc;
    speedI = TotalSpeedI / nt;  % average speed per steps
    
    % calculate the speed of the algorithm (over steps).


    resultI = {tip_trajI, whole_trajI, contact_trajI, param.d, speedI};

    results_mat{it} = resultI;

end

% the frictionless 

tip_traj_fl = zeros(3,nt);
whole_traj_fl = zeros(3*size(p,2), nt);

tip_traj_fl(:,1) = p(:,end);
whole_traj_fl(:,1) = reshape(p,[],1);

param_fl.d = 6;

tic

for i = 2:nt

    contacts_fl = [];
    [contacts_fl, ~, ~, ~] = detectAdditionalContacts(tube, u_fl, obstacles, contacts_fl, CornerRange1);
    % contact_traj{i-1} = contacts_fl;

    % actuation along X positive
    alphaI = alpha_traj(i);
    tube.T_base(1:2,4) = [x_traj(i), y_traj(i)]';  

    [u_fl, contacts_fl, ~, ~, p_fl] = getFrictionalContactShape3(tube, obstacles_fl, u_fl, contacts_fl, CornerRange2, param_fl);

    tip_traj(:,i) = p(:,end);
    whole_traj(:,i) = reshape(p,[],1);

    tip_traj_fl(:,i) = p_fl(:,end);
    whole_traj_fl(:,i) = reshape(p_fl,[],1);

end

TotalSpeedFl = toc;
speed_fl  = TotalSpeedFl / nt;
disp('Speed of frictionless model is ' +string(speed_fl*1000) +'ms.')

%% extract results

tip_traj_m6 = results_mat{1}{1};
tip_traj_m12 = results_mat{2}{1};
tip_traj_m150 = results_mat{end}{1};

whole_traj_m6 = results_mat{1}{2};
whole_traj_m12 = results_mat{2}{2};
whole_traj_m150 = results_mat{end}{2};

%% extract the force/contact traj

% contact1_force_traj = [];
% contact1_FricRatio_traj = [];
% for it = 1:length(contact_traj)
%     contactI = contact_traj(it);
%     contactI = contactI{1};
%     if length(contactI) == 1
%         fc = contactI(1).force;  % contact 1 is corner contact
%         nc = contactI(1).normal;
%         friction_ratio = norm((eye(3)-nc*nc')*fc) / norm(nc'*fc);
%         contact1_force_traj = [contact1_force_traj, fc];
%         contact1_FricRatio_traj = [contact1_FricRatio_traj, friction_ratio];
%     else
%         continue
%     end
% end
%%

obstacles2plot = obstacles{1};

figure('units','pixels','position',[0 0 700 400])
plotConfig3D([600 400])
hold on;
axis equal
grid on
box on

xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

view([2,30])
axis([-55, 200, -65, 65, -40, wall_dis+15])

pathCut = round(nt/3):nt;

plot3Mat(tip_traj_m150(:,pathCut), 'r', 1)
plot3Mat(tip_traj_m12(:,pathCut), [0.4660 0.6740 0.1880], 1)
plot3Mat(tip_traj_m6(:,pathCut), 'b', 1)
plot3Mat(tip_traj_fl(:,pathCut), 'k', 1)

% plot the actuation circle
h_act_traj = plot3Mat(act_traj(:,pathCut), 'k');
set(h_act_traj, 'LineStyle', '--');

%plot a representative shape
p_show = reshape(whole_traj_m150(:,round(nt/2)), 3, []);
plot3Mat(p_show, 'r', 2)

% plot the plane
plot_r = 65; 
[Xc,Yc,Zc] = cylinder(plot_r,100);
Zc2 = Zc*5 + wall_dis;
Xc2 = Xc + 140;
% surf(Xc2, Yc, Zc2)
h_env1 = surf(Xc2,Yc,Zc2, 'FaceColor', 'k',  'FaceAlpha', 0.2, 'EdgeColor', 0.4*[1 1 1], 'MeshStyle','row', 'LineWidth',0.5);
patchButt = fill3(Xc2(1,:),Yc(1,:),Zc2(1,:),'k','FaceAlpha', 0.2, 'EdgeColor','none');
patchTop = fill3(Xc2(2,:),Yc(2,:),Zc2(2,:),'k','FaceAlpha', 0.2, 'EdgeColor','none');

% saveas(gcf, "simu_outs/figs_case1/fig_case1_circle3D", 'png')
% close(gcf)

%% make vid

figure('units','pixels','position',[0 0 700 400])
plotConfig3D([600 400])
hold on;
axis equal
grid on
box on

xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')

view([2,45])
axis([-55, 200, -65, 65, -40, wall_dis+15])

h_m150 = plot3(nan,nan,nan, 'r', 'lineWidth', 2);
h_m12 = plot3(nan,nan,nan, 'Color', [0.4660 0.6740 0.1880], 'lineWidth', 2);
h_m6 = plot3(nan,nan,nan, 'b', 'lineWidth', 2);
h_fl = plot3(nan,nan,nan, 'k', 'lineWidth', 2);

htip_m150 = plot3(nan,nan,nan, 'r', 'lineWidth', 1);
htip_m12 = plot3(nan,nan,nan, 'Color', [0.4660 0.6740 0.1880], 'lineWidth', 1);
htip_m6 = plot3(nan,nan,nan, 'b', 'lineWidth', 1);
htip_fl = plot3(nan,nan,nan, 'k', 'lineWidth', 1);

h_base = plot3(nan,nan,nan, 'k--', 'LineWidth',1);

% plot environemnt
h_env1 = surf(Xc2,Yc,Zc2, 'FaceColor', 'k',  'FaceAlpha', 0.2, 'EdgeColor', 0.4*[1 1 1], 'MeshStyle','row', 'LineWidth',0.5);
patchButt = fill3(Xc2(1,:),Yc(1,:),Zc2(1,:),'k','FaceAlpha', 0.2, 'EdgeColor','none');
patchTop = fill3(Xc2(2,:),Yc(2,:),Zc2(2,:),'k','FaceAlpha', 0.2, 'EdgeColor','none');

pathCut = round(nt/3):nt;

n_s = length(pathCut);
%make vid
totalTime = 10;
vid3D_path = "simu_outs/figs_case1/vid_case1_circle";
vid = VideoWriter(vid3D_path,'Uncompressed AVI');
vid.FrameRate = round(n_s/totalTime);
open(vid)
i_start = pathCut(1);
for i = pathCut

    % draw shape

    pi_m150 = reshape(whole_traj_m150(:,i),3,[]);
    pi_m12 = reshape(whole_traj_m12(:,i),3,[]);
    pi_m6 = reshape(whole_traj_m6(:,i),3,[]);
    pi_fl = reshape(whole_traj_fl(:,i),3,[]);    
        
    set(h_m150, 'XData', pi_m150(1,:),'YData', pi_m150(2,:),'ZData', pi_m150(3,:))
    set(h_m12, 'XData', pi_m12(1,:),'YData', pi_m12(2,:),'ZData', pi_m12(3,:))
    set(h_m6, 'XData', pi_m6(1,:),'YData', pi_m6(2,:),'ZData', pi_m6(3,:))
    set(h_fl, 'XData', pi_fl(1,:),'YData', pi_fl(2,:),'ZData', pi_fl(3,:))


    % draw tip traj
    set(htip_m150, 'XData', tip_traj_m150(1,i_start:i),...
        'YData', tip_traj_m150(2,i_start:i),...
        'ZData', tip_traj_m150(3,i_start:i))
    set(htip_m12, 'XData', tip_traj_m12(1,i_start:i),...
        'YData', tip_traj_m12(2,i_start:i),...
        'ZData', tip_traj_m12(3,i_start:i))
    set(htip_m6, 'XData', tip_traj_m6(1,i_start:i),...
        'YData', tip_traj_m6(2,i_start:i),...
        'ZData', tip_traj_m6(3,i_start:i))
    set(htip_fl, 'XData', tip_traj_fl(1,i_start:i),...
        'YData', tip_traj_fl(2,i_start:i),...
        'ZData', tip_traj_fl(3,i_start:i))

    % draw base
    set(h_base, 'XData', act_traj(1,i_start:i),...
        'YData', act_traj(2,i_start:i),...
        'ZData', act_traj(3,i_start:i))

    drawnow
    writeVideo(vid,getframe(gcf));

end

close(vid)

%% calculation speed and accuracy.
% extract the speed
speed_list = [];
accuracy_list = [];
for i = 1:length(results_mat)
    tip_trajI = results_mat{i}{1};

    ni = length(pathCut);

    % compare
    [~, ~, ErrorsI] = compareShapes(tip_trajI(:,pathCut), tip_traj_m150(:,pathCut), 1:ni, false);

    speed_list = [speed_list, 1000*results_mat{i}{5}];

    accuracy_list = [accuracy_list, ErrorsI.maxBody];
end

figure()
plotConfig2D([300, 300])

subplot(2,1,1)
plotConfig2D([300, 300])
plot(d_list, speed_list, '-r.', 'LineWidth', 1)
ylabel('Time (ms)')

subplot(2,1,2)
plotConfig2D([300, 300])
plot(d_list, accuracy_list, '-r.', 'LineWidth', 1)
ylabel('Tip Error (mm)')
xlabel('Number of Edges')
ylim([0,max(accuracy_list)*1.05])

% saveas(gcf, "simu_outs/figs_case1/fig_case1_circle_speed", 'png')
% close(gcf)
%% 2D traj
figure()
plotConfig3D([300 300])
hold on
axis equal

plot(tip_traj_m150(1,pathCut), tip_traj_m150(2,pathCut), 'Color','r', 'LineWidth',1)
plot(tip_traj_m12(1,pathCut), tip_traj_m12(2,pathCut), 'Color',[0.4660 0.6740 0.1880], 'LineWidth',1)
plot(tip_traj_m6(1,pathCut), tip_traj_m6(2,pathCut), 'b', 'LineWidth',1)
plot(tip_traj_fl(1,pathCut), tip_traj_fl(2,pathCut), 'k', 'LineWidth',1)

xlim([85,195])
ylim([-55, 55])

% saveas(gcf, "simu_outs/figs_case1/fig_case1_circle_traj", 'png')
% close(gcf)

makeVid2D(gca, "simu_outs/figs_case1/vid_case1_circle_traj")
%% legend
figure()
plotConfig2D([600 200])
hold on
plot(nan, nan, 'Color','r', 'LineWidth',1)
plot(nan, nan, 'Color',[0.4660 0.6740 0.1880], 'LineWidth',1)
plot(nan, nan, 'b', 'LineWidth',1)
plot(nan, nan, 'k', 'LineWidth',1)
grid off

legend('m=150', 'm=12', 'm=6', 'Frictionless', 'NumColumns',4)
saveas(gcf, "simu_outs/figs_case1/fig_case1_circle_legend", 'png')
close(gcf)


