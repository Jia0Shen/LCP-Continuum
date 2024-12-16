function [u, contacts, T, R, p] = getInitialShapeEnergy(tube, obstacles, cornerFlag)
% Solve shape for tube-obstacle contacts
% constacts: array of struct(contact point, plane index, tube contact point index)
% obstacles: cells of obstacles

    function [E, dE] = potential_energy(u,uhat,du,K_mat)
        % delete the first element because we are doing numerical integral
        % K_mat(1:3,1:3) = zeros(3);
        if nargout > 1
            dE = du*K_mat*(u-uhat);
        end
        E = 0.5*du*(u-uhat)'*K_mat*(u-uhat);
    end

    function [c, ceq] = phi_constrain_u(tube, obstacles, u)
        u_3xn = reshape(u, 3, []);
        [T1, R1, p1] = solveShape(tube.T_base, u_3xn, tube.s);
        ceq = 0;
        c  =[];
        for io = 1:length(obstacles)
            ci = obstacles{io}.computeDepth(p1, false);
            % ci(isnan(ci)) = 0;
            % cih = - p1(3,:) + obstacles{io}.h;
            cih = obstacles{io}.computeCutHeight(p1);
            ci = min(ci, cih);
            c = [c; ci];
        end

    end

if nargin == 2
    cornerFlag = true;
end

    bPlot = false;

    nObs = length(obstacles);
    ns = length(tube.s);
    
    K = getTubeK(tube);
    K_mat = diag(K(:));
    % invK = 1./K;
    
    % init variables
    f = zeros(nObs, ns);       % plane contact forces
    u_hat = reshape(tube.uhat, [], 1);
    ds_step = (tube.s(end)-tube.s(1)) / (length(tube.s)-1);

    %use fmincon to solve for the minimum energy.

    options = optimoptions(@fmincon,...
                        'Algorithm', 'sqp', ...
                        'CheckGradients', false, ...
                        'SpecifyConstraintGradient', false, ...
                        'SpecifyObjectiveGradient', true, ...
                        "StepTolerance", 1e-6, ...
                        "ConstraintTolerance", 1e-10, ...
                        "OptimalityTolerance", 1e-5, ...
                        "MaxIteration", 2, ...
                        "MaxFunctionEvaluations", 1e10, ...
                        "OutputFcn", []);

                            % "EnableFeasibilityMode",true, ...
                        % "SubproblemAlgorithm","cg", ...

    A = [];lb = [];ub = [];b = [];Aeq = [];beq = [];
    obj_fun = @(u) potential_energy(u, u_hat, ds_step, K_mat);
    constraints = @(u) phi_constrain_u(tube, obstacles, u);
    u0_ini = u_hat; % zeros(size(u_hat));

    %solve and save
    [u_sol,Energy_min,~,~,Lagrangian] = fmincon(obj_fun,u0_ini,A,b,Aeq,beq,lb,ub,constraints,options); 

    [c_test, ~] = phi_constrain_u(tube, obstacles, u_sol);
    
    u = reshape(u_sol, 3, []);
    % [T, R, p] = solveShape(tube.T_base, u, tube.s);

    % for iter = 1:1
    %     contacts = [];
    %     [contacts, T, R, p] = detectAdditionalContacts(tube, u, obstacles, contacts, cornerFlag);
    %     [u, contacts, T, R, p] = getContactShape3(tube, obstacles, u, contacts, cornerFlag);
    % end

    contacts = [];
    [contacts, T, R, p] = detectAdditionalContacts(tube, u, obstacles, contacts, cornerFlag);

    % T, R, p
    if nargout > 2
        [T, R, p] = solveShape(tube.T_base, u, tube.s);
    end

end