function [u, contacts, T, R, p] = getContactShape3(tube, obstacles, u_prev, contacts_prev, cornerRange, maxFricFlag)
% Solve shape for tube-obstacle contacts
% constacts: array of struct(contact point, plane index, tube contact point index)
% obstacles: cells of obstacles

    if nargin == 5
        maxFricFlag = false;
    end

    [contacts_prev, T, R, p] = detectAdditionalContacts(tube, u_prev, obstacles, contacts_prev, cornerRange);

    if isempty(contacts_prev)
        contacts = contacts_prev;
        u = tube.uhat;
        [T, R, p] = solveShape(tube.T_base, u, tube.s);

        return;        
    end

    %% Find dp = Af + q
    % shape to linearize about
    [~, R_l, p_l] = solveShape(tube.T_base, u_prev, tube.s);
    
    J = computeJacobian(R_l, p_l);
    
    nc = length(contacts_prev);
    itc = zeros(1,nc);  % tube point ids on contact
    p_contacts = zeros(3, nc);
    obs_contacts_drift = zeros(3, nc);
    for ic = 1:nc
        itc(ic) = contacts_prev(ic).tube_point_id;
        p_contacts(:, ic) = contacts_prev(ic).tube_point;  %changed 3/20/2024, point or tube_point?
        obs_i = obstacles{contacts_prev(ic).obstacle_id};

        if size(obs_i.T_history, 3) >= 2
            T_curr = obs_i.T_history(:,:,end);
            T_prev = obs_i.T_history(:,:,end-1);
            pc_prev_i = contacts_prev(ic).point;
            pc_i = T_curr * inv(T_prev) * [pc_prev_i; 1];
            obs_contacts_drift(:, ic) = pc_i(1:3) - pc_prev_i;
        end
    end

    % Contact Jacobian
    itc3 = [3*itc-2; 3*itc-1; 3*itc];
    itc3 = itc3(:);
    
    Jc = J(itc3, :);

    % contact point drift by base transformation
    b = p_l(:,itc) - p_contacts - obs_contacts_drift;
    b = b(:);
    

    % prev moment
    K = getTubeK(tube);
    m_prev = K .* (u_prev(:) - tube.uhat(:));

    % ds
    % ds = tube.s(2) - tube.s(1);   % assumed constant
    ds = tube.s(end) - tube.s(end-1);   % assumed constant
    % A, q
    invK = 1./K;
    Jc_invK_ds = Jc .* invK' * ds;
    
    A = Jc_invK_ds * Jc';
    q = -Jc_invK_ds * m_prev + b;



    %% Construct LCP
    % n and mu
    n = zeros(3, nc, nc);
    mu = zeros(nc, nc);

    B = zeros(3,2,nc);  % contact tangent bases

    for ic = 1:nc
        io = contacts_prev(ic).obstacle_id;

        % n_cur = obstacles{io}.getNormal(p_contacts(:,ic));
        n_cur = contacts_prev(ic).normal;
        n(:,ic,ic) = n_cur;
        mu(ic,ic) = obstacles{io}.mu;

        B(:,:,ic) = null(n_cur');
    end
    n = reshape(n, [], nc);

    % D and e
    d = 6;          % number of edges in friction cone
    
    theta = linspace(0, pi, round(d/2)+1);
    theta = theta(1:end-1);

    cs = [cos(theta); sin(theta)];

    D = zeros(3*nc, d*nc);
    e = zeros(d*nc, nc);
    for ic = 1:nc
        D_cur = B(:,:,ic)*cs;
        D_cur = [D_cur, -D_cur];
        D(3*ic-2:3*ic, d*(ic-1)+1:d*ic) = D_cur;

        e(d*(ic-1)+1:d*ic, ic) = 1;
    end

    % keep all contacts
    ic_f = find(diag(mu) > -1e16);  % frictional contact ids
    ic_f = reshape(ic_f, 1, []);
    nc_f = length(ic_f);        % # of frictional contacts

    if nc_f ~= nc
        % mu, D, e
        mu = mu(ic_f, :);

        % row = 3*repmat(ic_f,3,1) - (2:-1:0)';
        col = d*repmat(ic_f,d,1) - (d-1:-1:0)';
        
        D = D(:, col);
        e = e(col, ic_f);
    end

    % construct LCP without friction terms

    % M = [n'*A*n, n'*A*D,   zeros(nc,nc_f);
    %      D'*A*n, D'*A*D,                e;
    %      mu,     -e',    zeros(nc_f,nc_f)];
    % 
    % g = [n'*q; D'*q; zeros(nc_f,1)];

    M = n'*A*n;
    g = n'*q;

    %% Solve LCP
    % [w,x] = LCPSolve(M,g);
    % 
    % fn = x(1);
    % beta = x(2:1+d);
    % lambda = x(end);

    % %% Solve scaled LCP
    sv = svd(A);
    scale = max(sv);

    dimM = size(M,1);

    Dw = diag([ones(1, dimM-nc_f), ones(1,nc_f)*scale]);
    Dx = diag([ones(1, dimM-nc_f)/scale, ones(1,nc_f)]);

    M2 = Dw * M * Dx;
    g2 = Dw * g;

    [w2, x2, retcode2] = LCPSolve(M2, g2, 1e-8, dimM^2);
    test_w2 = M2*x2+g;
    % retcode2  % 1 success; 2 max termination

    x = Dx * x2;
    % w = Dw \ w2;

    fn = x(1:nc);
    % beta = x(nc+1:nc+nc_f*d);
    % lambda = x(nc+nc_f*d+1:end);

    % disp([max(beta)/fn, lambda]);
    

    %% Solve for shape
    f = n*fn;
    % if ~isempty(beta)
    %     f = f + D*beta;
    % end
    e3 = [0 0 1]';
    f_fric = zeros(3,nc);
    if maxFricFlag
        % find the maximum friction direction
        % should have: tube.dalpha, tube.dbeta
        for ni = 1:nc
            id = contacts_prev(ni).tube_point_id;
            n_tangent = normalize(p(:,id) - p(:,id-1), 'norm');
            v_trans = tube.dbeta * n_tangent;
            v_rot = tube.dalpha * hat(e3) * contacts_prev(ni).tube_point;
            vt = v_trans + v_rot;
            dir_v = normalize(vt, 'norm');
            f_fric(:, ni) = - mu(ni,ni) * fn(ni) * dir_v;
        end
        
        % add friction.
        f_fric = reshape(f_fric, 3*nc, 1);
        f = f+ f_fric;
    end  

    m = Jc' * f; 
    u = reshape(invK .* m, 3, []) + tube.uhat;

    % T, R, p
    if nargout > 2
        [T, R, p] = solveShape(tube.T_base, u, tube.s);
    end

    %% Contacts
    contacts = contacts_prev;
    
    dp = A*f + q;
    p_contacts_new = p_contacts + reshape(dp, 3, []);

    ic_to_remove = [];

    for ic = 1:nc
        if (fn(ic) > 0)
            io = contacts(ic).obstacle_id;
            contacts(ic).point = obstacles{io}.project(p_contacts_new(:,ic));

            contacts(ic).force = f(3*ic-2:3*ic);
        else
            ic_to_remove = [ic_to_remove, ic];
        end
    end
    contacts(ic_to_remove) = [];

end