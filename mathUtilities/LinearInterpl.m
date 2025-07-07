function [y, x_idx] = LinearInterpl(x, dx)

% given a list x, interplate x linearly with step length dx.

n = length(x);
y = x(1);
idx = 1;

for i = 1:n-1

    if abs(x(i+1) - x(i)) < 2*dx
        addy = x(i+1);
        idx = [idx, 1];
        y = [y, addy];
    else

        ni = ceil(abs(x(i+1) - x(i)) / dx);

        addy = linspace(x(i), x(i+1), ni);

        addy(1) = [];

        y = [y, addy];
        idx = [idx, ni-1];
    end
end

x_idx = cumsum(idx);


end