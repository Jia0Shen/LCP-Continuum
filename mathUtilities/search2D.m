function [xsol, feval] = search2D(fun, lb, ub, x10, maxIter)
    if nargin == 4
        maxIter = 5;
    end
    options = optimset('MaxIter', 100);
    x1_sol = x10;
    
    for iter = 1:maxIter
        [x2_sol, feval] = fminbnd(@(x2) fun([x1_sol, x2]), lb(2), ub(2), options);
        [x1_sol, feval] = fminbnd(@(x1) fun([x1, x2_sol]), lb(1), ub(1), options);
    end

    xsol = [x1_sol, x2_sol];

end