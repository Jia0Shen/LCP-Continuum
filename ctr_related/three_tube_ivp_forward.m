function [g_end, shape, b, g_arrays] = three_tube_ivp_forward(ctr, q, w, yu)
    % forward kinematics with (q, w) returns the shape of the tube.
    % q = [alpha1, beta1, alpha2, beta2, alpha3, beta3];
    % w = [F, L];
    alpha1 = q(1); beta1 = q(2); alpha2 = q(3);
    beta2 = q(4); alpha3 = q(5); beta3 = q(6);
    F = w(1:3);  L = w(4:6);
    Ls3=ctr.Ls3; Ls2=ctr.Ls2; Ls1=ctr.Ls1;
    Lc3=ctr.Lc3; Lc2=ctr.Lc2; Lc1=ctr.Lc1;
    s1 = Ls3 + Lc3 + beta3;  % first section, 3 tubes
    s2 = Ls2 + Lc2 + beta2;  % second sectin, 2 tubes
    s3 = Ls1 + Lc1 + beta1;  % third section, 1 tube.
    E = ctr.E; G = ctr.G;
    ka1 = ctr.ka1; ka2 = ctr.ka2; ka3 = ctr.ka3;
    I1=ctr.I1; I2=ctr.I2; I3=ctr.I3; J1=ctr.J1; J2=ctr.J2; J3=ctr.J3;
    ode_solver = ctr.ode_solver;

    % yu = [u1z0, u2z0, u3z0, mbx0, mby0]
    u1z0 = yu(1); u2z0 = yu(2); u3z0 = yu(3);
    p0 = [0; 0; 0];
    R1_0 = [cos(alpha1-beta1*u1z0), -sin(alpha1-beta1*u1z0), 0;
            sin(alpha1-beta1*u1z0), cos(alpha1-beta1*u1z0),  0;
            0,                      0,                       1];
    theta2_0 = alpha2 - alpha1 - (beta2 * u2z0 - beta1 * u1z0);
    theta3_0 = alpha3 - alpha1 - (beta3 * u3z0 - beta1 * u1z0);
    ini_0_s1 = [p0; reshape(R1_0, 9, 1); 
                theta2_0; theta3_0; reshape(yu, 5,1)];
    [s_s1, y_s1] = ode_solver(@ode_0_s1, [0, s1], ini_0_s1);
    % y_s1, 3 tubes; y = [p,R,theta2,theta3,u1z,u2z,u3z,mbx,mby]  
    % y_s2, 2 tubes, y = [p,R,theta2,u1z,u2z,mbx,mby]
    % y_s3, 1 tube, y = [p,R,u1z,mbx,mby]
    ini_s1_s2 = y_s1(end, [1:13, 15, 16, 18, 19]);

    [s_s2, y_s2] = ode_solver(@ode_s1_s2, [s1,s2], ini_s1_s2);
    ini_s2_s3 = y_s2(end, [1:12, 14, 16, 17]);
    [s_s3, y_s3] = ode_solver(@ode_s2_s3, [s2, s3], ini_s2_s3);

    e3_main = [0;0;1];
    R1_main = reshape(y_s3(end, 4:12), 3, 3);
    u1z_l1 = y_s3(end, 13);
    u2z_l2 = y_s2(end, 15);
    u3z_l3 = y_s1(end, 17);
    mbxy_l1 = y_s3(end, 14:15)';
    b = [G * J1 * u1z_l1 - e3_main' * R1_main' * L;
         G * J2 * u2z_l2; 
         G * J3 * u3z_l3;
         mbxy_l1 - select_xy(R1_main' * L)];
    
    y_end = y_s3(end, :);
    state_var_end = y_end(13:15);
    p_end = y_end(1:3)';
    R_end = reshape(y_end(4:12), 3,3);
    g_end = [R_end, p_end; zeros(1,3), 1];

    % whole shape
    g_array1 = [y_s1(:,4:6),zeros(size(y_s1,1),1),...
        y_s1(:,7:9),zeros(size(y_s1,1),1),...
        y_s1(:,10:12),zeros(size(y_s1,1),1),...
        y_s1(:,1:3),ones(size(y_s1,1),1)];
    g_array2 = [y_s2(:,4:6),zeros(size(y_s2,1),1),...
        y_s2(:,7:9),zeros(size(y_s2,1),1),...
        y_s2(:,10:12),zeros(size(y_s2,1),1),...
        y_s2(:,1:3),ones(size(y_s2,1),1)];
    g_array3 = [y_s3(:,4:6),zeros(size(y_s3,1),1),...
        y_s3(:,7:9),zeros(size(y_s3,1),1),...
        y_s3(:,10:12),zeros(size(y_s3,1),1),...
        y_s3(:,1:3),ones(size(y_s3,1),1)];
    g_arrays = {g_array1, g_array2, g_array3};
    shape = {y_s1(:,1:3), y_s2(:,1:3), y_s3(:,1:3)};


    % if nargin == 4
    %     g_in = g_end;
    % else
    %     if s2 < s_in && s_in <= s3
    %         y_in = interp1(s_s3, y_s3, s_in);
    %         p_in = y_in(1:3)';
    %         R_in = reshape(y_in(4:12), 3,3);
    %         g_in = [R_in, p_in; zeros(1,3), 1];
    %     elseif s1 < s_in && s_in <= s2
    %         y_in = interp1(s_s2, y_s2, s_in);
    %         p_in = y_in(1:3)';
    %         R_in = reshape(y_in(4:12), 3,3);
    %         g_in = [R_in, p_in; zeros(1,3), 1];
    %     elseif 0 <= s_in && s_in <= s1
    %         y_in = interp1(s_s1, y_s1, s_in);
    %         p_in = y_in(1:3)';
    %         R_in = reshape(y_in(4:12), 3,3);
    %         g_in = [R_in, p_in; zeros(1,3), 1];
    %     else
    %         disp('no matched s_in')
    %     end
    % end

    
    function [ka1_s, ka2_s, ka3_s] = ka_s(s)
       % find the curvature at s (s>=0)
       % tube 3
       if s <= Ls3 + beta3
           ka3_s = ka3;  %  0
       elseif (s >= Ls3 + beta3) && (s <= s1)
           ka3_s = ka3;
       else
           ka3_s = nan;
       end

       % tube 2           
       if s <= Ls2 + beta2
           ka2_s = ka2;  %  0
       elseif (s >= Ls2 + beta2) && (s <= s2)
           ka2_s = ka2;
       else
           ka2_s = nan;
       end

       % tube 1
       if s <= Ls1 + beta1
           ka1_s = ka1;  %  0
       elseif (s >= Ls1 + beta1) && (s <= s3)
           ka1_s = ka1;
       else
           ka1_s = nan;
       end        
    end

    function dy = ode_0_s1(s, y)
        % 3 tubes; y = [p,R,theta2,theta3,u1z,u2z,u3z,mbx,mby]  
        [ka1_s, ka2_s, ka3_s] = ka_s(s);
        R1 = reshape(y(4:12), 3, 3);
        theta1 = 0;
        theta2 = y(13);  theta3 = y(14);
        u1z = y(15);  u2z = y(16);  u3z = y(17);
        mbx = y(18);  mby = y(19);

        mbz = G * (J1 * u1z + J2 * u2z + J3 * u3z);
        mb = [mbx; mby; mbz];
        sum_EIK = E * (I1*ka1_s*[cos(theta1);sin(theta1)] + ...
                       I2*ka2_s*[cos(theta2);sin(theta2)] + ...
                       I3*ka3_s*[cos(theta3);sin(theta3)] );
        u1xy = ([mbx; mby] + sum_EIK) / (E * (I1 + I2 + I3));
        u1y = u1xy(2);
        u2y = [- sin(theta2), cos(theta2)] * u1xy;
        u3y = [- sin(theta3), cos(theta3)] * u1xy;
        u1 = [u1xy; u1z];
        e3 = [0;0;1];

        dp = R1 * e3;
        dR1 = reshape(R1 * hat(u1), 9, 1);
        dtheta2 = u2z - u1z;
        dtheta3 = u3z - u1z;
        du1z = - (ka1_s * E * I1 / (G * J1)) * u1y;
        du2z = - (ka2_s * E * I2 / (G * J2)) * u2y;
        du3z = - (ka3_s * E * I3 / (G * J3)) * u3y;
        dmbxy = select_xy(- hat(u1) * mb - hat(e3) * R1' * F);

        dy = [dp; dR1; dtheta2; dtheta3; du1z; du2z; du3z; dmbxy]; 
    end

    function dy = ode_s1_s2(s, y)
        % 2 tubes, y = [p,R,theta2,u1z,u2z,mbx,mby]
        [ka1_s, ka2_s, ka3_s] = ka_s(s);
        R1 = reshape(y(4:12), 3, 3);
        theta1 = 0;
        theta2 = y(13);  
        u1z = y(14);  u2z = y(15);  
        mbx = y(16);  mby = y(17);

        mbz = G * (J1 * u1z + J2 * u2z);
        mb = [mbx; mby; mbz];
        sum_EIK = E * (I1*ka1_s*[cos(theta1);sin(theta1)] + ...
                       I2*ka2_s*[cos(theta2);sin(theta2)] );
        u1xy = ([mbx; mby] + sum_EIK) / (E * (I1 + I2));
        u1y = u1xy(2);
        u2y = [- sin(theta2), cos(theta2)] * u1xy;
        u1 = [u1xy; u1z];
        e3 = [0;0;1];

        dp = R1 * e3;
        dR1 = reshape(R1 * hat(u1), 9, 1);
        dtheta2 = u2z - u1z;
        du1z = - (ka1_s * E * I1 / (G * J1)) * u1y;
        du2z = - (ka2_s * E * I2 / (G * J2)) * u2y;
        dmbxy = select_xy(- hat(u1) * mb - hat(e3) * R1' * F);

        dy = [dp; dR1; dtheta2; du1z; du2z; dmbxy]; 

    end

    function dy = ode_s2_s3(s, y)
        % 1 tube y = [p,R,u1z,mbx,mby]
        [ka1_s, ka2_s, ka3_s] = ka_s(s);
        R1 = reshape(y(4:12), 3, 3);
        theta1 = 0;
        u1z = y(13);
        mbx = y(14);  mby = y(15);

        mbz = G * (J1 * u1z);
        mb = [mbx; mby; mbz];
        sum_EIK = E * (I1*ka1_s*[cos(theta1);sin(theta1)]);
        u1xy = ([mbx; mby] + sum_EIK) / (E * I1);
        u1y = u1xy(2);
        u1 = [u1xy; u1z];
        e3 = [0;0;1];

        dp = R1 * e3;
        dR1 = reshape(R1 * hat(u1), 9, 1);
        du1z = - (ka1_s * E * I1 / (G * J1)) * u1y;
        dmbxy = select_xy(- hat(u1) * mb - hat(e3) * R1' * F);

        dy = [dp; dR1; du1z; dmbxy]; 
    end
     
end 
