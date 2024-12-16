function h_env = plotConfigCylinder(width_height, obstacles, face_trans)

    function obstacles = buildEnv(PipeHeight, r, CoeffFrics)
        % obs1 is top, obs2 is bottom.
        z = [0;0;1];
        h1 = 1;
        h2 = 110 - h1;
        p1_cy = [0;0;PipeHeight-h1/2];
        p2_cy = [0;0;p1_cy(3)-(h1+h2)/2];
        obstacles{1} = objCylinderCut(p1_cy, z, r, h1, CoeffFrics(1), 'in');
        obstacles{2} = objCylinderCut(p2_cy, z, r, h2, CoeffFrics(2), 'in');
        obstacles{2}.cornerFlag = false;
    end

if nargin < 1
    width_height = [840 840];
    obstacles = buildEnv(160, 3.34, [0.01,0.01]);
elseif nargin == 1
    obstacles = buildEnv(160, 3.34, [0.01,0.01]);
elseif nargin == 2
    face_trans = 0.3;
end

    % biuld a new obstacle combining the two.
    zc = obstacles{1}.z;
    rc = obstacles{1}.r;
    p0c = obstacles{1}.p0 - 1/2*obstacles{2}.h*zc;
    hc = obstacles{1}.h + obstacles{2}.h;
    
    obstacle_c = objCylinderCut(p0c, zc, rc, hc, 0.01, 'in');

    [X,Y,Z] = obstacle_c.getMesh(200);
    h_env = surf(X,Y,Z, 'FaceColor', 'k',  'FaceAlpha', face_trans, 'EdgeColor', 'k', 'MeshStyle','row', 'EdgeAlpha',0.6, 'LineWidth',0.5);

    xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)'); 
    axis([-20, 50, -15, 15, 0, 200])

    axis equal
    grid on
    hold on

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