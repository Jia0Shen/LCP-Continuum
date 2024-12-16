function [J_ivp_one,C_ivp_one,E_end,V_ivp,B_end,J_in,C_in,g_end,b_end,E_in,g_in] = cal_Jacobian_ivp_one(ctr,q,w,yu0_true,s_in)
    % Calculate the jacobian using only one IVP.
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
        
    function dyVE = odeVE_0_s1(s, yVE)
        % 3 tubes; y = [p,R,theta2,theta3,u1z,u2z,u3z,mbx,mby,V,E]  
        [ka1_s, ka2_s, ka3_s] = ka_s(s);
        p = yVE(1:3);
        R1 = reshape(yVE(4:12), 3, 3);
        theta1 = 0;
        theta2 = yVE(13);  theta3 = yVE(14);
        u1z = yVE(15);  u2z = yVE(16);  u3z = yVE(17);
        mbx = yVE(18);  mby = yVE(19);
        % V = d/dxk[theta2,theta3,u1z,u2z,u3z,mbx,mby]
        VV = reshape(yVE(20:(20+7*17-1)), 7, 17);
        EE = reshape(yVE((20+7*17):(20+7*17+6*17-1)), 6, 17);
        
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

        dF_dx = [zeros(1,6), 1, 0, 0, zeros(1, 8);
                 zeros(1,6), 0, 1, 0, zeros(1, 8);
                 zeros(1,6), 0, 0, 1, zeros(1, 8)];
        g = [R1, p; zeros(1, 3), 1];
        % V = d/dxk[theta2,theta3,u1z,u2z,u3z,mbx,mby]
        
       sum_dEIK = E * (zeros(2, 17) + ...
           I2*ka2_s*[-sin(theta2)*ones(1,17);
                     cos(theta2)*ones(1,17)].*[VV(1,:);VV(1,:)] + ...
           I3*ka3_s*[-sin(theta3)*ones(1,17)
                     cos(theta3)*ones(1,17)].*[VV(2,:);VV(2,:)] );
       du1xy_dxk = 1 / (E * (I1 + I2 + I3)) * (VV(6:7, :) + sum_dEIK);
       du1_dxk = [du1xy_dxk; VV(3, :)];
       du1y_dxk = du1xy_dxk(2, :);
       du2y_dxk = VV(1,:) * ([-cos(theta2), - sin(theta2)] * u1xy) + ...
                    [-sin(theta2), cos(theta2)]*du1xy_dxk;
       du3y_dxk = VV(2,:) * ([-cos(theta3), - sin(theta3)] * u1xy) + ...
                    [-sin(theta3), cos(theta3)]*du1xy_dxk;

       % dR1_dxk = hat([zeros(3), eye(3)] * EE(:, k)) * R1;
       e3_dR1dxk_F = - hat(e3)*R1'*cross([zeros(3),eye(3)]*EE,F.*ones(3,17));
       dmb_dxk = [VV(6:7, :); G*(J1*VV(3,:)+J2*VV(4,:)+J3*VV(5,:))];
       ddmbxy_dxk = [eye(2), zeros(2,1)] * ...
                     (- cross(du1_dxk, mb.*ones(3,17)) ...
                      - cross(u1.*ones(3,17), dmb_dxk)  ...
                      - e3_dR1dxk_F ...
                      - hat(e3)*R1'*dF_dx);

       % creat dV
       dV = [
           VV(4,:) - VV(3,:);
           VV(5,:) - VV(3,:);
           - (ka1_s * E * I1) / (G * J1) * du1y_dxk;
           - (ka2_s * E * I2) / (G * J2) * du2y_dxk;
           - (ka3_s * E * I3) / (G * J3) * du3y_dxk;
           ddmbxy_dxk];

       % creat dE
       dE = Adj(g) * [zeros(3,17); du1_dxk];

       dV = reshape(dV, 7*17, 1);
       dE = reshape(dE, 6*17, 1);

       dyVE = [dp; dR1; dtheta2; dtheta3; ...
               du1z; du2z; du3z; dmbxy; dV; dE]; 
    end

    function dyVE = odeVE_s1_s2(s, yVE)
        % 2 tubes, y = [p,R,theta2,u1z,u2z,mbx,mby,V,E]
        [ka1_s, ka2_s, ka3_s] = ka_s(s);
        p = reshape(yVE(1:3),3,1);
        R1 = reshape(yVE(4:12), 3, 3);
        theta1 = 0;
        theta2 = yVE(13);  
        u1z = yVE(14);  u2z = yVE(15);  
        mbx = yVE(16);  mby = yVE(17);
        VV = reshape(yVE(18:(18+5*17-1)), 5, 17);
        EE = reshape(yVE((18+5*17):(18+5*17+6*17-1)), 6, 17);

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
        
        dV = zeros(5*17, 1);  dE = zeros(6*17,1);
        dF_dx = [zeros(1,6), 1, 0, 0, zeros(1, 8);
                 zeros(1,6), 0, 1, 0, zeros(1, 8);
                 zeros(1,6), 0, 0, 1, zeros(1, 8)];
        g = [R1, p; zeros(1, 3), 1];
        % V = d/dxk[theta2, u1z, u2z, mbx, mby];
        
       sum_dEIK = E * I2*ka2_s*[-sin(theta2)*ones(1,17);
                 cos(theta2)*ones(1,17)].*[VV(1,:);VV(1,:)];
       du1xy_dxk = 1 / (E * (I1 + I2)) * (VV(4:5, :) + sum_dEIK);
       du1_dxk = [du1xy_dxk; VV(2, :)];
       du1y_dxk = du1xy_dxk(2, :);
       du2y_dxk = VV(1,:) * ([-cos(theta2), - sin(theta2)]*u1xy) + ...
                            [-sin(theta2), cos(theta2)]*du1xy_dxk;
       % dR1_dxk = hat([zeros(3), eye(3)] * EE(:, k)) * R1;
       e3_dR1dxk_F = - hat(e3)*R1'*cross([zeros(3),eye(3)]*EE,F.*ones(3,17));
       dmb_dxk = [VV(4:5, :); G*(J1*VV(2,:)+J2*VV(3,:))];
       ddmbxy_dxk = [eye(2), zeros(2,1)] * ...
                     (- cross(du1_dxk, mb.*ones(3,17)) ...
                      - cross(u1.*ones(3,17), dmb_dxk)  ...
                      - e3_dR1dxk_F ...
                      - hat(e3)*R1'*dF_dx);

       % creat dV
       dV = [
           VV(3,:) - VV(2,:);
           - (ka1_s * E * I1) / (G * J1) * du1y_dxk;
           - (ka2_s * E * I2) / (G * J2) * du2y_dxk;
           ddmbxy_dxk];

       % creat dE
       dE = Adj(g) * [zeros(3,17); du1_dxk];
        
        dV = reshape(dV, 5*17, 1);
        dE = reshape(dE, 6*17, 1);
        
        dyVE = [dp; dR1; dtheta2; du1z; du2z; dmbxy; dV; dE]; 

    end

    function dyVE = odeVE_s2_s3(s, yVE)
        % 1 tube y = [p,R,u1z,mbx,mby,V,E]
        [ka1_s, ka2_s, ka3_s] = ka_s(s);
        p = reshape(yVE(1:3), 3,1);
        R1 = reshape(yVE(4:12), 3, 3);
        theta1 = 0;
        u1z = yVE(13);
        mbx = yVE(14);  mby = yVE(15);
        VV = reshape(yVE(16:(16+3*17-1)), 3, 17);
        EE = reshape(yVE((16+3*17):(16+3*17+6*17-1)), 6, 17);

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
        
        dV = zeros(3*17, 1);  dE = zeros(6*17,1);
        dF_dx = [zeros(1,6), 1, 0, 0, zeros(1, 8);
                 zeros(1,6), 0, 1, 0, zeros(1, 8);
                 zeros(1,6), 0, 0, 1, zeros(1, 8)];
        g = [R1, p; zeros(1, 3), 1];
        % V = d/dxk[u1z, mbx, mby]
        
        % sum_dEIK = E * (I1*ka1_s*[-sin(theta1);cos(theta1)] * 0);
        sum_dEIK = 0;
        du1xy_dxk = 1 / (E * I1) * (VV(2:3, :) + sum_dEIK);
        du1_dxk = [du1xy_dxk; VV(1, :)];
        du1y_dxk = du1xy_dxk(2, :);
        % dR1_dxk = hat([zeros(3), eye(3)] * EE(:, k)) * R1;
        % dF_dxk = dF_dx(:, k);
        dmb_dxk = [VV(2:3, :); G * J1 * VV(1, :)];
        e3_dR1dxk_F = - hat(e3)*R1'*cross([zeros(3),eye(3)]*EE,F.*ones(3,17));
        
        ddmbxy_dxk = [eye(2), zeros(2,1)] * ...
                     (- cross(du1_dxk, mb.*ones(3,17)) ...
                      - cross(u1.*ones(3,17), dmb_dxk)  ...
                      - e3_dR1dxk_F ...
                      - hat(e3)*R1'*dF_dx);
        % creat dV
        dV = [
           - (ka1_s * E * I1) / (G * J1) * du1y_dxk;
           ddmbxy_dxk];

        % creat dE
        dE = Adj(g) * [zeros(3,17); du1_dxk];
        
        dV = reshape(dV, 3*17, 1);
        dE = reshape(dE, 6*17, 1);

        dyVE = [dp; dR1; du1z; dmbxy; dV; dE]; 
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

    function [V_end, E_end, B_end, V_in, E_in, g_end, b_end, g_in] = int_VE(qq, ww, yu, s_in)
        alpha1 = qq(1); beta1 = qq(2); alpha2 = qq(3);
        beta2 = qq(4); alpha3 = qq(5); beta3 = qq(6);
        F = ww(1:3);  L = ww(4:6);
        s1 = Ls3 + Lc3 + beta3;  % first section, 3 tubes
        s2 = Ls2 + Lc2 + beta2;  % second sectin, 2 tubes
        s3 = Ls1 + Lc1 + beta1;  % third section, 1 tube.
        % yu = [u1z0, u2z0, u3z0, mbx0, mby0], get V0, E0 from yu
        u1z0 = yu(1); u2z0 = yu(2); u3z0 = yu(3);
        p0 = [0; 0; 0];
        R1_0 = [cos(alpha1-beta1*u1z0), -sin(alpha1-beta1*u1z0), 0;
                sin(alpha1-beta1*u1z0), cos(alpha1-beta1*u1z0),  0;
                0,                      0,                       1];
        theta2_0 = alpha2 - alpha1 - (beta2 * u2z0 - beta1 * u1z0);
        theta3_0 = alpha3 - alpha1 - (beta3 * u3z0 - beta1 * u1z0);
        g0 = [R1_0, p0; zeros(1,3), 1];
        
        V0 = [-1,u1z0,1,-u2z0,0,    0,zeros(1,6),beta1,-beta2,     0,0,0;
              -1,u1z0,0,    0,1,-u3z0,zeros(1,6),beta1,     0,-beta3,0,0;
              zeros(5,6),             zeros(5,6),                 eye(5)];
        E0_c1 = [0,0,0,0,0,1]';
        E0_c2 = [0,0,0,0,0,-u1z0]';
        E0_c13 = [0,0,0,0,0,-beta1]';
        E0 = [E0_c1, E0_c2, zeros(6,10), E0_c13, zeros(6,4)];
        
        % 0 -> s1
        ini_0_s1 = [p0; reshape(R1_0, 9, 1); 
                    theta2_0; theta3_0; reshape(yu, 5,1);
                    reshape(V0, 7*17,1); reshape(E0, 6*17,1)];
        [s_s1, y_s1] = ode_solver(@odeVE_0_s1, [0, s1], ini_0_s1);
        
        % y_s1, 3 tubes; y = [p,R,theta2,theta3,u1z,u2z,u3z,mbx,mby]  
        % y_s2, 2 tubes, y = [p,R,theta2,u1z,u2z,mbx,mby]
        % y_s3, 1 tube, y = [p,R,u1z,mbx,mby]
        % ini_s1_s2 = y_s1(end, [1:13, 15, 16, 18, 19]);
        V_s1 = reshape(y_s1(end, 20:(20+7*17-1)), 7, 17);

        E_s1 = reshape(y_s1(end, (20+7*17):(20+7*17+6*17-1)), 6, 17);
        g_s1 = [reshape(y_s1(end, 4:12), 3,3), y_s1(end, 1:3)'; zeros(1,3),1];
        dy_s1_all = ode_0_s1(s1, y_s1(end, 1:19));
        dy_s1_ = dy_s1_all(13:19);
        dg_s1_ = [reshape(dy_s1_all(4:12),3,3), dy_s1_all(1:3); zeros(1,3), 0];
        
        dy_s1_plus_all = ode_s1_s2(s1, y_s1(end, [1:13, 15, 16, 18, 19]));
        dy_s1_plus = dy_s1_plus_all(13:17);
        dg_s1_plus = [reshape(dy_s1_plus_all(4:12),3,3), dy_s1_plus_all(1:3);zeros(1,3),0];
        
        % y_s1, 3 tubes; y = [theta2,theta3,u1z,u2z,u3z,mbx,mby]  
        % y_s2, 2 tubes, y = [theta2,u1z,u2z,mbx,mby]
        % y_s3, 1 tube, y = [u1z,mbx,mby]
        V_s1_plus = V_s1([1,3,4,6,7],:);
        V_s1_plus(:,6) = V_s1_plus(:,6) ...
                           + dy_s1_([1,3,4,6,7],:) - dy_s1_plus;
        V_s1_full_der = V_s1;
        V_s1_full_der(:, 6) = V_s1_full_der(:,6) + dy_s1_ ;
%                                 - [dy_s1_plus(1);
%                                    0;
%                                    dy_s1_plus(2:3);
%                                    0; 
%                                    dy_s1_plus(4:5)];
        E_s1_plus = E_s1;
        E_s1_plus(:, 6) = E_s1_plus(:, 6) ...
                                + inv_hat_se(dg_s1_ / g_s1) ...
                                - inv_hat_se(dg_s1_plus / g_s1);
        E_s1_full_der = E_s1;
        E_s1_full_der(:, 6) = E_s1_full_der(:, 6) ...
                                + inv_hat_se(dg_s1_ / g_s1);
        % % b / A = b * inv(A);
        % aa_test = Adj(g_s1) * inv_hat_se(inv(g_s1) * dg_s1_) - inv_hat_se(dg_s1_ * inv(g_s1));
        
        % s1 -> s2
        % ini_s2_s3 = y_s2(end, [1:12, 14, 16, 17]);
        ini_s1_s2 = [y_s1(end, [1:13, 15, 16, 18, 19]), ...
                        reshape(V_s1_plus, 1, 5*17), ...
                        reshape(E_s1_plus, 1, 6*17)];
        [s_s2, y_s2] = ode_solver(@odeVE_s1_s2, [s1,s2], ini_s1_s2);
        
        V_s2 = reshape(y_s2(end, 18:(18+5*17-1)), 5, 17);
        E_s2 = reshape(y_s2(end, (18+5*17):(18+5*17+6*17-1)), 6, 17);
        
        g_s2 = [reshape(y_s2(end, 4:12), 3,3), y_s2(end, 1:3)'; zeros(1,3),1];
        dy_s2_all = ode_s1_s2(s2, y_s2(end, 1:17));
        dy_s2_ = dy_s2_all(13:17);
        dg_s2_ = [reshape(dy_s2_all(4:12),3,3), dy_s2_all(1:3); zeros(1,3), 0];
        
        dy_s2_plus_all = ode_s2_s3(s2, y_s2(end, [1:12,14,16,17]));
        dy_s2_plus = dy_s2_plus_all(13:15);
        dg_s2_plus = [reshape(dy_s2_plus_all(4:12),3,3), dy_s2_plus_all(1:3);zeros(1,3),0];
        
        % y_s1, 3 tubes; y = [theta2,theta3,u1z,u2z,u3z,mbx,mby]  
        % y_s2, 2 tubes, y = [theta2,u1z,u2z,mbx,mby]
        % y_s3, 1 tube, y = [u1z,mbx,mby]
        V_s2_plus = V_s2([2,4,5],:);
        V_s2_plus(:,4) = V_s2_plus(:,4) + ...
                            dy_s2_([2,4,5]) - dy_s2_plus;
        V_s2_full_der = V_s2;
        V_s2_full_der(:, 4) = V_s2_full_der(:,4) + dy_s2_;
%                                 - [0;
%                                    dy_s2_plus(1);
%                                    0;
%                                    dy_s2_plus(2:3)];

        E_s2_plus = E_s2;                                   
        E_s2_plus(:, 4) = E_s2_plus(:, 4) ...
                                + inv_hat_se(dg_s2_ / g_s2) ...
                                - inv_hat_se(dg_s2_plus / g_s2);
        E_s2_full_der = E_s2;                                   
        E_s2_full_der(:, 4) = E_s2_full_der(:, 4) ...
                                + inv_hat_se(dg_s2_ / g_s2);
        ini_s2_s3 = [y_s2(end, [1:12, 14, 16, 17]), ...
                        reshape(V_s2_plus, 1, 3*17), ...
                        reshape(E_s2_plus, 1, 6*17)];
        [s_s3, y_s3] = ode_solver(@odeVE_s2_s3, [s2, s3], ini_s2_s3);
        V_s3 = reshape(y_s3(end, 16:(16+3*17-1)), 3, 17);
        E_s3 = reshape(y_s3(end, (16+3*17):(16+3*17+6*17-1)), 6, 17);
        R1_l1 = reshape(y_s3(end, 4:12), 3, 3);

        % s3 -> end transition
        g_s3 = [reshape(y_s3(end, 4:12), 3,3), y_s3(end, 1:3)'; zeros(1,3),1];
        dy_s3_all = ode_s2_s3(s3, y_s3(end, 1:15));
        dy_s3_ = dy_s3_all(13:15);
        dg_s3_ = [reshape(dy_s3_all(4:12),3,3), dy_s3_all(1:3); zeros(1,3), 0];
        V_s3_plus = V_s3;
        E_s3_plus = E_s3;
        % aa_test = Adj(g_s1) * inv_hat_se(inv(g_s1) * dg_s1_) - inv_hat_se(dg_s1_ * inv(g_s1));
        V_s3_plus(:,2) = V_s3_plus(:,2) + dy_s3_;
        E_s3_plus(:,2) = E_s3_plus(:,2) ...
                                + inv_hat_se(dg_s3_ / g_s3);
        
        V_end = V_s3_plus;
        E_spatial = E_s3_plus;
%         V = V_s3;
%         E = E_s3;
        
        p_end = g_s3(1:3, 4);
        % the value we interest in
        if norm(s_in - s3) < 1e-15
            V_in = V_end;
            E_in_spatial = E_spatial;
            p_in = p_end;
            g_in = g_s3;
        elseif s2 < s_in && s_in < s3
            y_in = interp1(s_s3, y_s3, s_in);
            p_in = y_in(1:3);
            R_in = reshape(y_in(4:12), 3,3);
            g_in = [R_in, p_in'; zeros(1,3), 1];
            V_in = reshape(y_in(end, 16:(16+3*17-1)), 3, 17);
            E_in_spatial = reshape(y_in(end, (16+3*17):(16+3*17+6*17-1)), 6, 17);
        elseif s1 < s_in && s_in <= s2
            y_in = interp1(s_s2, y_s2, s_in);
            p_in = y_in(1:3);
            R_in = reshape(y_in(4:12), 3,3);
            g_in = [R_in, p_in'; zeros(1,3), 1];
            V_in = reshape(y_in(end, 18:(18+5*17-1)), 5, 17);
            E_in_spatial = reshape(y_in(end, (18+5*17):(18+5*17+6*17-1)), 6, 17);
        elseif 0 <= s_in && s_in <= s1
            y_in = interp1(s_s1, y_s1, s_in);
            p_in = y_in(1:3);
            R_in = reshape(y_in(4:12), 3,3);
            g_in = [R_in, p_in'; zeros(1,3), 1];
            V_in = reshape(y_in(:, 20:(20+7*17-1)), 7, 17);
            E_in_spatial = reshape(y_in(:, (20+7*17):(20+7*17+6*17-1)), 6, 17);
        else
            V_in = V_end;
            E_in_spatial = E_spatial;
            p_in = p_end;
            g_in = g_s3;
            disp('no matched s_in')
        end
        

        % convert to hybrid Jacobian
        p_end = g_s3(1:3, 4);
        E_end = [eye(3), -hat(p_end); zeros(3,3), eye(3)] * E_spatial;
        E_in = [eye(3), -hat(p_in); zeros(3,3), eye(3)] * E_in_spatial;

        % y_s1, 3 tubes; y = [theta2,theta3,u1z,u2z,u3z,mbx,mby]  
        % y_s2, 2 tubes, y = [theta2,u1z,u2z,mbx,mby]
        % y_s3, 1 tube, y = [u1z,mbx,mby]
        % calculate B from V, E
        % e3_dR1dxk_F = - hat(e3)*R1'*cross([zeros(3),eye(3)]*EE,F.*ones(3,17));
        dL_dxk = [zeros(3, 9), eye(3), zeros(3, 5)];
        RT_dL_dxk = R1_l1'*dL_dxk;
        dRT_dkx_L = - R1_l1' * cross([zeros(3),eye(3)]*E_end, L.*ones(3,17));
        e3 = [0,0,1]';

        u1z_l1 = y_s3(end, 13);
        u2z_l2 = y_s2(end, 15);
        u3z_l3 = y_s1(end, 17);
        mbxy_l1 = y_s3(end, 14:15)';

        g_end = g_s3;

        b_end = [G * J1 * u1z_l1 - e3' * R1_l1' * L;
             G * J2 * u2z_l2; 
             G * J3 * u3z_l3;
             mbxy_l1 - select_xy(R1_l1' * L)];

        B_end = [
            G*J1*V_end(1,:) - e3'*(dRT_dkx_L + RT_dL_dxk);
            G*J2*V_s2_full_der(3, :);
            G*J3*V_s1_full_der(5, :);
            V_end(2:3, :) - dRT_dkx_L(1:2, :) - RT_dL_dxk(1:2, :);
        ];
    end     

    if nargin == 4
        s_in = ctr.L1 + q(2);
    end
    [V_ivp, E_end, B_end, V_in, E_in, g_end, b_end, g_in] = int_VE(q, w, yu0_true, s_in);

    % get jacobian from E, B.
    Eu_end = E_end(:, 13:17);
    Eq_end = E_end(:, 1:6);  
    Ew_end = E_end(:, 7:12);
    Bu_end = B_end(:, 13:17);
    Bq_end = B_end(:, 1:6);
    Bw_end = B_end(:, 7:12);

    Eu_in = E_in(:, 13:17);
    Eq_in = E_in(:, 1:6);  
    Ew_in = E_in(:, 7:12);

    
    J_ivp_one = Eq_end - Eu_end * inv(Bu_end) * Bq_end;
    C_ivp_one = Ew_end - Eu_end * inv(Bu_end) * Bw_end;

    J_in = Eq_in - Eu_in * inv(Bu_end) * Bq_end;
    C_in = Ew_in - Eu_in * inv(Bu_end) * Bw_end;
    
end
