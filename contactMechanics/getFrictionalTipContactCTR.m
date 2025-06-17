function [state,contactFlag] = getFrictionalTipContactCTR(ctr, obs, q, state0, d)
% solve the updated shape if the tip is in contact with obstacle.
    if nargin == 4
        d = 6;    % tunable
    end
    p_prev = state0.p;
    w_prev = state0.w;
    yu0_prev = state0.yu0;
    q_prev = state0.q;

    % using the previous Jacobian
    [J, C, ~, ~, B] = cal_Jacobian_ivp_one(ctr, q_prev, w_prev, yu0_prev);

    Jc = J(1:3,:);
    Cc = C(1:3, 1:3);   % Cc = dp/df

    % geometry
    p_tip = p_prev{3}(end,:)';
    n = obs.getNormal(p_tip);
    mu = obs.mu;

    % d = 6;   % tunable
    theta = linspace(0, pi, round(d/2)+1);
    theta = theta(1:end-1);
    cs = [cos(theta); sin(theta)];
    nullVec = null(n');
    D = [nullVec*cs, -nullVec*cs];
    e = ones(d,1);

    % fixed the penetration/vibration
    p_proj = obs.project(p_tip);
    b = - (p_proj - p_tip);

    f_prev = w_prev(1:3);
    h = Jc*(q-q_prev) - Cc*f_prev+ b;

    % LCP
    M = [n'*Cc*n, n'*Cc*D, 0;
         D'*Cc*n, D'*Cc*D, e;
         mu,      -e',     0];

    g = [n'*h; D'*h; 0];
    
    % solve scaled LCP
    sv = svd(Cc);
    scale = max(sv);
    dimM = size(M,1);

    Dw = diag([ones(1, dimM-1), scale]);
    Dx = diag([ones(1, dimM-1)/scale, 1]);

    M2 = Dw * M * Dx;
    g2 = Dw * g;

    [w2, x2, retcode2] = LCPSolve(M2, g2, 1e-8, dimM^2);
    test_w2 = M2*x2+g2;
    x_sol = Dx * x2;

    fn = x_sol(1);
    beta = x_sol(2:d+1);

    % calculate the update state
    f = n*fn + D*beta;
    w = [f; zeros(3,1)];

    % solve the bvp using new w
    % [~, p, yu0, b_res, shape] = three_tube_fk(ctr, q, w, yu0_prev);
    % state.p = p;
    % state.w = w;
    % state.yu0 = yu0;
    % state.q = q;
    % state.g = shape;

    % solve the yu0 using Jacobian update
    Bu = B(:, 13:17);
    Bq = B(:, 1:6);
    Bw = B(:, 7:12);

    dyu0 = - Bu \ (Bq*(q-q_prev) + Bw*(w-w_prev));  % see Caleb's paper
    yu0 = yu0_prev + dyu0;

    % use this for faster speed
    % [g_end, p, b_res, shape] = three_tube_ivp_forward(ctr, q, w, yu0);

    % use this for better stability
    [~, p, yu0, b_res, shape] = three_tube_fk(ctr, q, w, yu0);

    state.p = p;
    state.w = w;
    state.yu0 = yu0;
    state.q = q;
    state.g = shape;

    % check the acc of compliance
    dp = Jc*(q-q_prev) + Cc*(f-f_prev);
    dp_num = p{3}(end,:)' - p_tip;
    % check if the contact is maintained.
    if fn > 1e-12
        contactFlag = true;
    else 
        contactFlag = false;
    end


end