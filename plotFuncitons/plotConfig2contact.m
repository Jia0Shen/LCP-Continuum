function plotConfig2contact(width_height, obstacles)

    function obstacles = buildEnv(h_corn, h_wall, mu_s)
        z = [0;0;1];
        width_corner = 5;
        p1_cy = [0;0;h_corn-width_corner/2];
        r = 5 - 0.66;  % tube outer radii is 0.66mm.
        
        obstacles{1} = objCylinderCut(p1_cy, z, r, width_corner, mu_s(1), 'in');
        obstacles{1}.cornerFlag = true;
    
        p0_plane = [0;0;h_wall];
        z_plane = [0;0;-1];
        
        obstacles{2} = objPlane(p0_plane, z_plane, mu_s(2));
    end

if nargin < 1
    width_height = [840 840];
    obstacles = buildEnv(125, 170, 0.01);
elseif nargin ==1
    obstacles = buildEnv(125, 170, 0.01);
end

    [X1,Y1,Z1] = obstacles{1}.getMesh(50);
    [X2,Y2,Z2] = obstacles{2}.getMesh(80, 100);
    surf(X1,Y1,Z1, 'FaceColor', 'k',  'FaceAlpha', 0.3, 'EdgeColor', 'k', 'MeshStyle','row', 'EdgeAlpha',0.6, 'LineWidth',0.5);
    surf(X2,Y2,Z2, 'FaceColor', 'k',  'FaceAlpha', 0.3, 'EdgeColor', 'k', 'MeshStyle','row', 'EdgeAlpha',0.6, 'LineWidth',0.5);
    fill3(X2(1,:),Y2(1,:),Z2(1,:),'k','FaceAlpha', 0.3);
    fill3(X2(2,:),Y2(2,:),Z2(2,:),'k','FaceAlpha', 0.3);

    xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); 
    axis([-20, 100, -85, 85, 0, 200])

    set(gcf, 'Units', 'pixels');
    set(gcf, 'position', [0 0 width_height]);
    % paper size
    set(gcf, 'PaperPositionMode', 'auto');
    %white background
    set(gcf, 'Color', 'w');

    t=get(gca,'Title');
    set(t,'FontSize',11);    % 18 
    
    %set axis label font size
    tx = get(gca, 'Xlabel');
    ty = get(gca, 'Ylabel');
    tz = get(gca, 'Zlabel');
    set(tx,'FontSize',11);    % 16
    set(ty,'FontSize',11);
    set(tz,'FontSize',11);
    
    % set axis text size and legend text size
    set(gca,'FontSize',11);   % 14

    set(gca, ...
        'Box'         , 'on'     , ...
        'TickDir'     , 'out'     , ...
        'TickLength'  , [.01 .01] , ...
        'XMinorTick'  , 'on'      , ...
        'YMinorTick'  , 'on'      , ...
        'XGrid'       , 'on'      , ...
        'YGrid'       , 'on'      , ...
        'LineWidth'   ,  1); 

end