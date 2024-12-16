function h = plot3Mat(p, color, width)

if nargin == 1
    h = plot3(p(1,:), p(2,:), p(3,:));
elseif nargin == 2
    h = plot3(p(1,:), p(2,:), p(3,:), 'Color', color);
elseif nargin == 3
    h = plot3(p(1,:), p(2,:), p(3,:), 'Color', color, 'LineWidth',width);
end
    
end