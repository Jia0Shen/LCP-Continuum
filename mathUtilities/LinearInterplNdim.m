function [y, x_idx] = LinearInterplNdim(x, dx)

% given a list x m by n: each column is a vector of point.

[m,n] = size(x);
y = x(:,1);
idx = 1;

for i = 1:n-1

    if norm(x(:,i+1) - x(:,i)) < 2*dx
        addy = x(:,i+1);
        idx = [idx, 1];
        y = [y, addy];
    else

        ni = ceil(norm(x(:,i+1) - x(:,i)) / dx);

        % addy = linspace(x(i), x(i+1), ni);
        addy = x(:,i) + (x(:,i+1) - x(:,i))*(1:ni-1)/(ni-1);

        % addy(1) = [];

        y = [y, addy];
        idx = [idx, ni-1];
    end
end

x_idx = cumsum(idx);


end