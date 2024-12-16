function plotConfig2D(width_height)

if nargin < 1
    width_height = [840 840];
end

    hold on
    grid on
    % axis equal
    % xlabel('x (mm)'); ylabel('y (mm)');

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
    set(tx,'FontSize',11);    % 16
    set(ty,'FontSize',11);
    
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