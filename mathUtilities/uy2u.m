function u = uy2u(uy)
    n = length(uy);
    uy = reshape(uy, 1, n);
    u = [zeros(1, n); uy; zeros(1, n)];
    u = reshape(u, 3*n, 1);
end