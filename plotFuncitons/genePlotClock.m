function [tick, circle] = genePlotClock(height, R, r, alpha, grids)

if nargin == 4
    grids = 12;
end

angles = linspace(0,2*pi,grids+1) + alpha;
angles_dense = linspace(0,2*pi,100) + alpha;

tick.x = [R*sin(angles); r*sin(angles)];
tick.y = [R*cos(angles); r*cos(angles)];
tick.z = ones(size(tick.x))*height;

circle.x = [R*sin(angles_dense)', r*sin(angles_dense)'];
circle.y = [R*cos(angles_dense)', r*cos(angles_dense)'];
circle.z = ones(size(circle.x))*height;


end