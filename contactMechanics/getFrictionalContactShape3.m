function [u, contacts, T, R, p] = getFrictionalContactShape3(tube, obstacles, u_prev, contacts_prev, cornerRange, param)
% Solve shape for tube-obstacle contacts
% constact: (contact point, contact force, obstacle index, tube contact point index)
% contact: + contact.type; contact.normal; contact.tube_point
% **** For now, only implemented for constant ds ****

if nargin == 5
    d = 150;  % number of edges in friction cone
elseif nargin > 5
    d = param.d;   % number of edges in friction cone
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
    % itc = zeros(1,nc);  % tube point ids on contact
    pc_vec = zeros(3, nc); % location of the projection of contact point after obstacle moves
    % obs_contacts_drift = zeros(3, nc);
    Jc = zeros(3*nc,size(J,2));
    % p_lc = zeros(3, nc); 
    b_vec = zeros(3, nc);

    for ic = 1:nc
        itc = contacts_prev(ic).tube_point_id;
        pt = contacts_prev(ic).tube_point;
        % p_contacts(:, ic) = pt;  %changed 3/20/2024, point or tube_point?
        obs_i = obstacles{contacts_prev(ic).obstacle_id};

        if strcmp(contacts_prev(ic).type, 'cornerContact')
            itc3_1 = [3*itc(1)-2; 3*itc(1)-1; 3*itc(1)];
            itc3 = [3*itc(2)-2; 3*itc(2)-1; 3*itc(2)];
            % find distance from p_(i-1), pc, and p_i
            pi_1 = pt(:,1);  p_i = pt(:,2);
            p_prev = contacts_prev(ic).point;
            di_1 = norm(pi_1-p_prev); di = norm(p_i-p_prev);
            ds = di+di_1;

            % linearly combine the jacobian
            Jc(3*ic-2:3*ic, :) = 1/ds*(di*J(itc3_1, :) + di_1*J(itc3,:));
            % p_contacts(:, ic) = pc;  % pc is already projected to the surface.
            % p_proj = obs_i.project2corner(p_prev, cornerRange/2);
            p_proj = obs_i.project2corner(p_prev, 0);
            p_lc = 1/ds*(di*p_l(:,itc(1)) + di_1*p_l(:,itc(2)));
        
        else
            itc3 = [3*itc-2; 3*itc-1; 3*itc];
            p_proj = obs_i.project(pt);
            p_lc = p_l(:, itc);
            Jc(3*ic-2:3*ic, :) = J(itc3, :);
        end

        if size(obs_i.T_history, 3) >= 2
            T_curr = obs_i.T_history(:,:,end);
            T_prev = obs_i.T_history(:,:,end-1);
            % pc_prev_i = contacts_prev(ic).point;
            pc_i = T_curr * inv(T_prev) * [p_proj; 1];
            pc_vec(:,ic) = pc_i(1:3);
            % obs_contacts_drift(:, ic) = pc_i(1:3) - pc_prev_i;
        else
            % no obstacle movement
            pc_vec(:,ic) = p_proj;
        end

        b_vec(:,ic) = p_lc - pc_vec(:,ic);
    end

    % contact point drift by base transformation
    b = b_vec(:);
    
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

    % keep frictional contacts only
    ic_f = find(diag(mu) > 0);  % frictional contact ids
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

    M = [n'*A*n, n'*A*D,   zeros(nc,nc_f);
         D'*A*n, D'*A*D,                e;
         mu,     -e',    zeros(nc_f,nc_f)];

    g = [n'*q; D'*q; zeros(nc_f,1)];

    %% Solve LCP
    % [w,x] = LCPSolve(M,g);
    % 
    % fn = x(1);
    % beta = x(2:1+d);
    % lambda = x(end);


    % %% Solve scaled LCP
    sv = svd(A);
    scale = 1; % max(sv);

    dimM = size(M,1);

    Dw = diag([ones(1, dimM-nc_f), ones(1,nc_f)*scale]);
    Dx = diag([ones(1, dimM-nc_f)/scale, ones(1,nc_f)]);

    M2 = Dw * M * Dx;
    g2 = Dw * g;

    [w2, x2, retcode2] = LCPSolve(M2, g2, 1e-8, dimM^2);
    test_w2 = M2*x2+g2;

    % disp('relative diff of w2/w2_ is ' + string(norm(w2-test_w2)/norm(w2)));
    % [w2, x2] = LCPSolve(M2, g2, 1e-16);
    % [w2, x2] = LCP_jh(M2, g2);

    x = Dx * x2;
    % w = Dw \ w2;

    fn = x(1:nc);
    beta = x(nc+1:nc+nc_f*d);
    % lambda = x(nc+nc_f*d+1:end);

    % disp([max(beta)/fn, lambda]);
    

    %% Solve for shape
    f = n*fn;
    if ~isempty(beta)
        f = f + D*beta;
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
    test = n'*dp;
    p_contacts_new = pc_vec + reshape(dp, 3, []);

    ic_to_remove = [];

    for ic = 1:nc
        if (fn(ic) > 0)
            io = contacts(ic).obstacle_id;
            contacts(ic).force = f(3*ic-2:3*ic);
            % contacts(ic).point = obstacles{io}.project(p_contacts_new(:,ic));
            contacts(ic).point = p_contacts_new(:,ic);

        else
            ic_to_remove = [ic_to_remove, ic];
        end
    end
    contacts(ic_to_remove) = [];

end













