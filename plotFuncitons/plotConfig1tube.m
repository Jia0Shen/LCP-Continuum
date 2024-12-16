function plotConfig1tube(width_height, obstacles)

    function obstacles = buildEnv(PlaneHeight, mu)
        p0_plane = [0;0;PlaneHeight];
        z_plane = [0;0;-1];
        obstacles{1} = objPlane(p0_plane, z_plane, mu);
    end

if nargin < 1
    width_height = [840 840];
    obstacles = buildEnv(160, 0.01);
elseif nargin ==1
    obstacles = buildEnv(160, 0.01);
end
    % 'units','pixels','position',[0 0 840 840]
    % draw obstacles
    [X1,Y1,Z1] = obstacles{1}.getMesh(100, 100);
    surf(X1,Y1,Z1, 'FaceColor', 'k',  'FaceAlpha', 0.3, 'EdgeColor', 'k', 'MeshStyle','row', 'EdgeAlpha',0.6, 'LineWidth',0.5);
    fill3(X1(1,:),Y1(1,:),Z1(1,:),'k','FaceAlpha', 0.3);
    fill3(X1(2,:),Y1(2,:),Z1(2,:),'k','FaceAlpha', 0.3);

    xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); 
    axis([-15, 110, -100, 100, 0, 180])

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