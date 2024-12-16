function [hess,jacob_compl,hess_in,jacob_compl_in,J_end,C_end,J_in,C_in,E_end,B_end,DB,D_end,g_end,b_end,pJ_px,pC_px] = cal_Hessian_ivpforward(ctr, q, w, yu0_true, s_in)
    % pJ_px,pC_px,U_end,D_end,DB,V_end,E_end,B_end
    function [ka1_s, ka2_s, ka3_s] = ka_s(s)
        ka1_s = ka1;
        ka2_s = ka2;
        ka3_s = ka3;
    end

    function dyH = ode_Hess_3_tube(s, yH)
        % 3 tubes; y = [p,R,theta2,theta3,u1z,u2z,u3z,mbx,mby,V,E,U,D]  
        % r = 1;
        [ka1_s, ka2_s, ka3_s] = ka_s(s);
        p = yH(1:3); R = reshape(yH(4:12),3,3);
        theta2 = yH(13);  theta3 = yH(14);
        u1z = yH(15);  u2z = yH(16);  u3z = yH(17);
        mbx = yH(18);  mby = yH(19);
        % V = d/dxk[theta2,theta3,u1z,u2z,u3z,mbx,mby]
        mbz = G * (J1 * u1z + J2 * u2z + J3 * u3z);
        mb = [mbx; mby; mbz];
        sum_EIK = E * (I1*ka1_s*[1;0] + ...
                       I2*ka2_s*[cos(theta2);sin(theta2)] + ...
                       I3*ka3_s*[cos(theta3);sin(theta3)] );
        u1xy = ([mbx; mby] + sum_EIK) / (E * (I1 + I2 + I3));
        u1y = u1xy(2);
        u2y = [- sin(theta2), cos(theta2)] * u1xy;
        u3y = [- sin(theta3), cos(theta3)] * u1xy;
        u1 = [u1xy; u1z];

        e3 = [0;0;1];
        dp = R * e3;
        dR = R * hat(u1); dR_vec = reshape(dR, 9, 1);
        dtheta2 = u2z - u1z;
        dtheta3 = u3z - u1z;
        du1z = - (ka1_s * E * I1 / (G * J1)) * u1y;
        du2z = - (ka2_s * E * I2 / (G * J2)) * u2y;
        du3z = - (ka3_s * E * I3 / (G * J3)) * u3y;
        dmbxy = select_xy(- hat(u1) * mb - hat(e3)*R'*F);

        if(length(yH) == 19)
            dyH = [dp; dR_vec; dtheta2; dtheta3; du1z; du2z; du3z; dmbxy];
            return
        end

        % V = d/dxk[theta2,theta3,u1z,u2z,u3z,mbx,mby]

        O3I3 = [zeros(3,3), eye(3)];
        VV = reshape(yH(20:(20+7*17-1)), 7, 17);
        EE = reshape(yH(20+7*17: 20+7*17+6*17-1), 6, 17);
        dF_dx = [zeros(3,6), eye(3), zeros(3,8)];
        
        sum_dEIK = E * (zeros(2, 17) + ...
           I2*ka2_s*[-sin(theta2)*ones(1,17);
                     cos(theta2)*ones(1,17)].*[VV(1,:);VV(1,:)] + ...
           I3*ka3_s*[-sin(theta3)*ones(1,17)
                     cos(theta3)*ones(1,17)].*[VV(2,:);VV(2,:)] );
        pu1xy_pxk = 1 / (E * (I1 + I2 + I3)) * (VV(6:7, :) + sum_dEIK);
        pu1_pxk = [pu1xy_pxk; VV(3, :)];
        pu1y_pxk = pu1xy_pxk(2, :);
        pu2y_pxk = VV(1,:) * ([-cos(theta2), - sin(theta2)] * u1xy) + ...
                    [-sin(theta2), cos(theta2)]*pu1xy_pxk;
        pu3y_pxk = VV(2,:) * ([-cos(theta3), - sin(theta3)] * u1xy) + ...
                    [-sin(theta3), cos(theta3)]*pu1xy_pxk;
        pmb_pxk = [VV(6:7, :); G*(J1*VV(3,:)+J2*VV(4,:)+J3*VV(5,:))];
        e3_dR1dxk_F = - hat(e3)*R'*cross([zeros(3),eye(3)]*EE,F.*ones(3,17));
        d_pmbxy_pxk_ds = [eye(2), zeros(2,1)] * ...
                     (- cross(pu1_pxk, mb.*ones(3,17)) ...
                      - cross(u1.*ones(3,17), pmb_pxk) ...
                      - e3_dR1dxk_F ...
                      - hat(e3)*R'*dF_dx);
        % creat dV dE
        dV = [
           VV(4,:) - VV(3,:);
           VV(5,:) - VV(3,:);
           - (ka1_s * E * I1) / (G * J1) * pu1y_pxk;
           - (ka2_s * E * I2) / (G * J2) * pu2y_pxk;
           - (ka3_s * E * I3) / (G * J3) * pu3y_pxk;
           d_pmbxy_pxk_ds];

        dE = [-hat(R*e3) * O3I3 * EE;
             R * pu1_pxk];
        
        dV_vec = reshape(dV, 7*17, 1);
        dE_vec = reshape(dE, 6*17, 1);

        if(length(yH)==19+6*17+7*17)
             dyH = [dp; dR_vec; dtheta2; dtheta3; ...
               du1z; du2z; du3z; dmbxy; dV_vec; dE_vec];
             return
        end

        % second derivative
        dU = zeros(7,17,17);  dD = zeros(6,17,17);
        UU = reshape(yH(20+7*17+6*17: 20+7*17+6*17+7*17*17-1), 7, 17, 17);
        DD = reshape(yH(20+7*17+6*17+7*17*17: 20+7*17+6*17+7*17*17+6*17*17-1),6,17,17);

        [~, kkk, rrr] = size(dU); 
        % rr = 6; % for shrinked state
        ppu1xy_pxkpxr = 1/(E*(I1+I2+I3)) * (UU(6:7,:,:) + E*( ...
           I2*ka2_s*pagemtimes([-sin(theta2);cos(theta2)],UU(1,:,:)) + ...
           I2*ka2_s*pagemtimes([-cos(theta2);-sin(theta2)], tensorlize(VV(1,:)'*VV(1,1:rrr))) + ...
           I3*ka3_s*pagemtimes([-sin(theta3);cos(theta3)], UU(2,:,:)) + ...
           I3*ka3_s*pagemtimes([-cos(theta3);-sin(theta3)], tensorlize(VV(2,:)'*VV(2,1:rrr))) ));
        ppu1_pxkpxr = [ppu1xy_pxkpxr;UU(3,:,:)];
        ppu2y_pxkpxr = ([-cos(theta2),-sin(theta2)]*u1xy)*UU(1,:,:) + ...
           ([sin(theta2),-cos(theta2)]*u1xy)*tensorlize(VV(1,:)'*VV(1,1:rrr)) + ...
           tensorlize(([-cos(theta2),-sin(theta2)]*pu1xy_pxk)'*VV(1,1:rrr)) + ...
           tensorlize((VV(1,:)'*([-cos(theta2),-sin(theta2)]*pu1xy_pxk(:,1:rrr)))) + ...
           pagemtimes([-sin(theta2),cos(theta2)], ppu1xy_pxkpxr);
        ppu3y_pxkpxr = ([-cos(theta3),-sin(theta3)]*u1xy)*UU(2,:,:) + ...
           ([sin(theta3),-cos(theta3)]*u1xy)*tensorlize(VV(2,:)'*VV(2,1:rrr)) + ...
           tensorlize(([-cos(theta3),-sin(theta3)]*pu1xy_pxk)'*VV(2,1:rrr)) + ...
           tensorlize((VV(2,:)'*([-cos(theta3),-sin(theta3)]*pu1xy_pxk(:,1:rrr)))) + ...
           pagemtimes([-sin(theta3),cos(theta3)], ppu1xy_pxkpxr);
        ppmbz_pxkpxr = G*(J1*UU(3,:,:) + J2*UU(4,:,:) + J3*UU(5,:,:));
        ppmb_pxkpxr = [UU(6:7,:,:); ppmbz_pxkpxr];
                
        d_ppmbxy_pxkpxr_ds = select_xy( ...
           pagemtimes(hat(mb),ppu1_pxkpxr) ...  % to here
           + cross(repmat(reshape(pmb_pxk(:,1:rrr),[3,1,rrr]),[1,kkk,1]) , repmat(pu1_pxk, [1,1,rrr])) ...
           - cross(repmat(reshape(pu1_pxk(:,1:rrr),[3,1,rrr]),[1,kkk,1]) , repmat(pmb_pxk, [1,1,rrr])) ... % hat(pu1_pxk(:,r))*pmb_pxk ...
           - pagemtimes(hat(u1), ppmb_pxkpxr) ...
           + pagemtimes(pagemtimes(hat(e3)*R', hat_mat(O3I3*EE(:,1:rrr))), hat(F)*O3I3*(EE)+dF_dx) ...
           - pagemtimes(hat(e3)*R'*hat(F)*O3I3, DD(:,:,1:rrr)) ...
           );
        
        % find dU, dD
    
        dU(:,:,:) = cat(1, ...
           UU(4,:,:) - UU(3,:,:), ...
           UU(5,:,:) - UU(3,:,:), ...
           - (ka1_s * E * I1) / (G * J1) * ppu1xy_pxkpxr(2,:,:), ...
           - (ka2_s * E * I2) / (G * J2) * ppu2y_pxkpxr, ...
           - (ka3_s * E * I3) / (G * J3) * ppu3y_pxkpxr, ...
           d_ppmbxy_pxkpxr_ds);
    
        dD(:,:,1:rrr) = cat(1, ...
            - pagemtimes(hat(R*e3)*O3I3, DD(:,:,1:rrr)) ...
            + pagemtimes(hat_mat(- pagemtimes(hat_mat(O3I3*EE(:,1:rrr)), R*e3)), O3I3*EE), ...
              pagemtimes(hat_mat(O3I3*EE(:,1:rrr)), R*pu1_pxk) ...
              + pagemtimes(R, ppu1_pxkpxr));

        dU_vec = reshape(dU, 7*17*17, 1);
        dD_vec = reshape(dD, 6*17*17, 1);
        
        dyH = [dp; dR_vec; dtheta2; dtheta3; ...
               du1z; du2z; du3z; dmbxy; dV_vec; dE_vec; dU_vec; dD_vec]; 
    end

    function dyH = ode_Hess_2_tube(s, yH)
        % 2 tubes; y = [p,R,theta2,u1z,u2z,mbx,mby,V,E,U,D]  
        % r = 1;
        [ka1_s, ka2_s, ka3_s] = ka_s(s);
        p = yH(1:3); R = reshape(yH(4:12),3,3);
        theta2 = yH(13);  u1z = yH(14); u2z = yH(15);  
        mbx = yH(16);  mby = yH(17);
        % V = d/dxk[theta2,u1z,u2z,mbx,mby]

        mbz = G * (J1 * u1z + J2 * u2z);
        mb = [mbx; mby; mbz];
        sum_EIK = E * (I1*ka1_s*[1;0] + ...
                       I2*ka2_s*[cos(theta2);sin(theta2)]);
        u1xy = ([mbx; mby] + sum_EIK) / (E * (I1 + I2));
        u1y = u1xy(2);
        u2y = [- sin(theta2), cos(theta2)] * u1xy;
        u1 = [u1xy; u1z];

        e3 = [0;0;1];
        dp = R * e3;
        dR = R * hat(u1); dR_vec = reshape(dR, 9, 1);
        dtheta2 = u2z - u1z;
        du1z = - (ka1_s * E * I1 / (G * J1)) * u1y;
        du2z = - (ka2_s * E * I2 / (G * J2)) * u2y;
        dmbxy = select_xy(- hat(u1) * mb - hat(e3)*R'*F);

        if(length(yH)==17)
            dyH = [dp; dR_vec; dtheta2; du1z; du2z; dmbxy];
            return
        end

        % V = d/dxk[theta2, u1z, u2z, mbx, mby];
        O3I3 = [zeros(3,3), eye(3)];
        VV = reshape(yH(18:(18+5*17-1)), 5, 17);
        EE = reshape(yH(18+5*17: 18+5*17+6*17-1), 6, 17);
        dF_dx = [zeros(3,6), eye(3), zeros(3,8)];

        sum_dEIK = E*I2*ka2_s*[-sin(theta2)*ones(1,17);
                     cos(theta2)*ones(1,17)].*[VV(1,:);VV(1,:)];
        pu1xy_pxk = 1 / (E * (I1 + I2)) * (VV(4:5, :) + sum_dEIK);
        pu1_pxk = [pu1xy_pxk; VV(2, :)];
        pu1y_pxk = pu1xy_pxk(2, :);
        pu2y_pxk = VV(1,:) * ([-cos(theta2), - sin(theta2)] * u1xy) + ...
                    [-sin(theta2), cos(theta2)]*pu1xy_pxk;
        pmb_pxk = [VV(4:5, :); G*(J1*VV(2,:)+J2*VV(3,:))];
        e3_dR1dxk_F = - hat(e3)*R'*cross([zeros(3),eye(3)]*EE,F.*ones(3,17));
        d_pmbxy_pxk_ds = [eye(2), zeros(2,1)] * ...
                     (- cross(pu1_pxk, mb.*ones(3,17)) ...
                      - cross(u1.*ones(3,17), pmb_pxk) ...
                      - e3_dR1dxk_F ...
                      - hat(e3)*R'*dF_dx);
        % creat dV dE
        dV = [
           VV(3,:) - VV(2,:);
           - (ka1_s * E * I1) / (G * J1) * pu1y_pxk;
           - (ka2_s * E * I2) / (G * J2) * pu2y_pxk;
           d_pmbxy_pxk_ds];

        dE = [-hat(R*e3) * O3I3 * EE;
             R * pu1_pxk];        
        
        dV_vec = reshape(dV, 5*17, 1);
        dE_vec = reshape(dE, 6*17, 1);

        if(length(yH)==17+5*17+6*17)
            dyH = [dp; dR_vec; dtheta2; du1z; du2z; dmbxy;dV_vec;dE_vec];
            return
        end

        % second derivative
        dU = zeros(5,17,17);  dD = zeros(6,17,17);
        UU = reshape(yH(18+5*17+6*17: 18+5*17+6*17+5*17*17-1), 5, 17, 17);
        DD = reshape(yH(18+5*17+6*17+5*17*17: 18+5*17+6*17+5*17*17+6*17*17-1),6,17,17);

        % ----------------- vectorize -------------------
        [~, kkk, rrr] = size(dU);
        % rr = 6; % for shrinked state
        ppu1xy_pxkpxr = 1/(E*(I1+I2)) * (UU(4:5,:,:) + E*( ...
           I2*ka2_s*pagemtimes([-sin(theta2);cos(theta2)],UU(1,:,:)) + ...
           I2*ka2_s*pagemtimes([-cos(theta2);-sin(theta2)], tensorlize(VV(1,:)'*VV(1,1:rrr))) ));
        ppu1_pxkpxr = [ppu1xy_pxkpxr;UU(2,:,:)];
        ppu2y_pxkpxr = ([-cos(theta2),-sin(theta2)]*u1xy)*UU(1,:,:) + ...
           ([sin(theta2),-cos(theta2)]*u1xy)*tensorlize(VV(1,:)'*VV(1,1:rrr)) + ...
           tensorlize(([-cos(theta2),-sin(theta2)]*pu1xy_pxk)'*VV(1,1:rrr)) + ...
           tensorlize((VV(1,:)'*([-cos(theta2),-sin(theta2)]*pu1xy_pxk(:,1:rrr)))) + ...
           pagemtimes([-sin(theta2),cos(theta2)], ppu1xy_pxkpxr);
        ppmbz_pxkpxr = G*(J1*UU(2,:,:) + J2*UU(3,:,:));
        ppmb_pxkpxr = [UU(4:5,:,:); ppmbz_pxkpxr];
    
        d_ppmbxy_pxkpxr_ds = select_xy( ...
           pagemtimes(hat(mb),ppu1_pxkpxr) ...  % to here
           + cross(repmat(reshape(pmb_pxk(:,1:rrr),[3,1,rrr]),[1,kkk,1]) , repmat(pu1_pxk, [1,1,rrr])) ...
           - cross(repmat(reshape(pu1_pxk(:,1:rrr),[3,1,rrr]),[1,kkk,1]) , repmat(pmb_pxk, [1,1,rrr])) ... % hat(pu1_pxk(:,r))*pmb_pxk ...
           - pagemtimes(hat(u1), ppmb_pxkpxr) ...
           + pagemtimes(pagemtimes(hat(e3)*R', hat_mat(O3I3*EE(:,1:rrr))), hat(F)*O3I3*(EE)+dF_dx) ...
           - pagemtimes(hat(e3)*R'*hat(F)*O3I3, DD(:,:,1:rrr)) ...
           );
        
        % find dU, dD
    
        dU(:,:,:) = cat(1, ...
           UU(3,:,:) - UU(2,:,:), ...
           - (ka1_s * E * I1) / (G * J1) * ppu1xy_pxkpxr(2,:,:), ...
           - (ka2_s * E * I2) / (G * J2) * ppu2y_pxkpxr, ...
           d_ppmbxy_pxkpxr_ds);
    
        dD(:,:,1:rrr) = cat(1, ...
            - pagemtimes(hat(R*e3)*O3I3, DD(:,:,1:rrr)) ...
            + pagemtimes(hat_mat(- pagemtimes(hat_mat(O3I3*EE(:,1:rrr)), R*e3)), O3I3*EE), ...
              pagemtimes(hat_mat(O3I3*EE(:,1:rrr)), R*pu1_pxk) ...
              + pagemtimes(R, ppu1_pxkpxr));

        dU_vec = reshape(dU, 5*17*17, 1);
        dD_vec = reshape(dD, 6*17*17, 1);
        
        dyH = [dp; dR_vec; dtheta2; du1z; du2z; dmbxy; 
               dV_vec; dE_vec; dU_vec; dD_vec]; 
    end

    function dyH = ode_Hess_1_tube(s, yH)
        % 1 tubes; y = [p,R,u1z,mbx,mby,V,E,U,D]  
        % r = 1;
        [ka1_s, ka2_s, ka3_s] = ka_s(s);
        p = yH(1:3); R = reshape(yH(4:12),3,3);
        u1z = yH(13);  
        mbx = yH(14);  mby = yH(15);
        % V = d/dxk[u1z,mbx,mby]

        mbz = G * J1 * u1z;
        mb = [mbx; mby; mbz];
        sum_EIK = E * (I1*ka1_s*[1;0]);
        u1xy = ([mbx; mby] + sum_EIK) / (E * I1);
        u1y = u1xy(2);
        u1 = [u1xy; u1z];

        e3 = [0;0;1];
        dp = R * e3;
        dR = R * hat(u1); dR_vec = reshape(dR, 9, 1);
        du1z = - (ka1_s * E * I1 / (G * J1)) * u1y;
        dmbxy = select_xy(- hat(u1) * mb - hat(e3)*R'*F);

        if (length(yH)==15)
            dyH = [dp; dR_vec; du1z; dmbxy];
            return
        end

        % V = d/dxk[theta2,theta3,u1z,u2z,u3z,mbx,mby]
        O3I3 = [zeros(3,3), eye(3)];
        VV = reshape(yH(16:(16+3*17-1)), 3, 17);
        EE = reshape(yH(16+3*17: 16+3*17+6*17-1), 6, 17);
        dF_dx = [zeros(3,6), eye(3), zeros(3,8)];

        pu1xy_pxk = 1 / (E * I1) * VV(2:3, :);
        pu1_pxk = [pu1xy_pxk; VV(1, :)];
        pu1y_pxk = pu1xy_pxk(2, :);
        pmb_pxk = [VV(2:3, :); G*J1*VV(1,:)];
        e3_dR1dxk_F = - hat(e3)*R'*cross([zeros(3),eye(3)]*EE,F.*ones(3,17));
        d_pmbxy_pxk_ds = [eye(2), zeros(2,1)] * ...
                     (- cross(pu1_pxk, mb.*ones(3,17)) ...
                      - cross(u1.*ones(3,17), pmb_pxk) ...
                      - e3_dR1dxk_F ...
                      - hat(e3)*R'*dF_dx);
        % creat dV dE
        dV = [
           - (ka1_s * E * I1) / (G * J1) * pu1y_pxk;
           d_pmbxy_pxk_ds];

        dE = [-hat(R*e3) * O3I3 * EE;
             R * pu1_pxk];        
        
        dV_vec = reshape(dV, 3*17, 1);
        dE_vec = reshape(dE, 6*17, 1);

        if (length(yH)==15+3*17+6*17)
            dyH = [dp; dR_vec; du1z; dmbxy; dV_vec;dE_vec];
            return
        end

        % second derivative
        dU = zeros(3,17,17);  dD = zeros(6,17,17);
        UU = reshape(yH(16+3*17+6*17: 16+3*17+6*17+3*17*17-1), 3, 17, 17);
        DD = reshape(yH(16+3*17+6*17+3*17*17: 16+3*17+6*17+3*17*17+6*17*17-1),6,17,17);

        % ----------------- vectorize -------------------
        [~, kkk, rrr] = size(dU);
        % rr = 6; % for shrinked state
        ppu1xy_pxkpxr = 1/(E*I1) * UU(2:3,:,:);
        ppu1_pxkpxr = [ppu1xy_pxkpxr;UU(1,:,:)];
        ppmbz_pxkpxr = G*(J1*UU(1,:,:));
        ppmb_pxkpxr = [UU(2:3,:,:); ppmbz_pxkpxr];
    
        d_ppmbxy_pxkpxr_ds = select_xy( ...
           pagemtimes(hat(mb),ppu1_pxkpxr) ...  % to here
           + cross(repmat(reshape(pmb_pxk(:,1:rrr),[3,1,rrr]),[1,kkk,1]) , repmat(pu1_pxk, [1,1,rrr])) ...
           - cross(repmat(reshape(pu1_pxk(:,1:rrr),[3,1,rrr]),[1,kkk,1]) , repmat(pmb_pxk, [1,1,rrr])) ... % hat(pu1_pxk(:,r))*pmb_pxk ...
           - pagemtimes(hat(u1), ppmb_pxkpxr) ...
           + pagemtimes(pagemtimes(hat(e3)*R', hat_mat(O3I3*EE(:,1:rrr))), hat(F)*O3I3*(EE)+dF_dx) ...
           - pagemtimes(hat(e3)*R'*hat(F)*O3I3, DD(:,:,1:rrr)) ...
           );
        
        % find dU, dD
    
        dU(:,:,:) = cat(1, ...
           - (ka1_s * E * I1) / (G * J1) * ppu1xy_pxkpxr(2,:,:), ...
           d_ppmbxy_pxkpxr_ds);
    
        dD(:,:,1:rrr) = cat(1, ...
            - pagemtimes(hat(R*e3)*O3I3, DD(:,:,1:rrr)) ...
            + pagemtimes(hat_mat(- pagemtimes(hat_mat(O3I3*EE(:,1:rrr)), R*e3)), O3I3*EE), ...
              pagemtimes(hat_mat(O3I3*EE(:,1:rrr)), R*pu1_pxk) ...
              + pagemtimes(R, ppu1_pxkpxr));

        dU_vec = reshape(dU, 3*17*17, 1);
        dD_vec = reshape(dD, 6*17*17, 1);
        
        dyH = [dp; dR_vec;  du1z; dmbxy; 
            dV_vec; dE_vec; dU_vec; dD_vec]; 
    end
        
    function [p,R,state_7,V,E,U,D] = unzip_3_tube(y_3997)
        p = y_3997(1:3); R = reshape(y_3997(4:12),3,3);
        state_7 = y_3997(13:19);
        V = nan; E = nan; U = nan; D = nan;
        if (length(y_3997)==19)
            return
        end

        % V = d/dxk[theta2,theta3,u1z,u2z,u3z,mbx,mby]
        V = reshape(y_3997(20:(20+7*17-1)), 7, 17);
        E = reshape(y_3997(20+7*17: 20+7*17+6*17-1), 6, 17);
        if (length(y_3997) == 19+7*17+6*17)
            return
        end

        U = reshape(y_3997(20+7*17+6*17: 20+7*17+6*17+7*17*17-1), 7, 17, 17);
        D = reshape(y_3997(20+7*17+6*17+7*17*17: 20+7*17+6*17+7*17*17+6*17*17-1),6,17,17);
    end

    function [p,R,state_5,V,E,U,D] = unzip_2_tube(y_3383)
        p = y_3383(1:3); R = reshape(y_3383(4:12),3,3);
        state_5 = y_3383(13:17);
        V = nan; E = nan; U = nan; D = nan;
        if (length(y_3383)==17)
            return
        end

        % V = d/dxk[theta2,theta3,u1z,u2z,u3z,mbx,mby]
        V = reshape(y_3383(18:(18+5*17-1)), 5, 17);
        E = reshape(y_3383(18+5*17: 18+5*17+6*17-1), 6, 17);

        if (length(y_3383)==17+5*17+6*17)
            return
        end

        U = reshape(y_3383(18+5*17+6*17: 18+5*17+6*17+5*17*17-1), 5, 17, 17);
        D = reshape(y_3383(18+5*17+6*17+5*17*17: 18+5*17+6*17+5*17*17+6*17*17-1),6,17,17);
    end

    function [p,R,state_3,V,E,U,D] = unzip_1_tube(y_2769)
        p = y_2769(1:3); R = reshape(y_2769(4:12),3,3);
        state_3 = y_2769(13:15);
        V = nan; E = nan; U = nan; D = nan;

        if (length(y_2769)==15)
            return
        end

        % V = d/dxk[theta2,theta3,u1z,u2z,u3z,mbx,mby]
        V = reshape(y_2769(16:(16+3*17-1)), 3, 17);
        E = reshape(y_2769(16+3*17: 16+3*17+6*17-1), 6, 17);

        if (length(y_2769)==15+3*17+6*17)
            return
        end

        U = reshape(y_2769(16+3*17+6*17: 16+3*17+6*17+3*17*17-1), 3, 17, 17);
        D = reshape(y_2769(16+3*17+6*17+3*17*17: 16+3*17+6*17+3*17*17+6*17*17-1),6,17,17);
    end

    function y_5 = state_725(y_7)
        if (length(size(y_7)) < 3)
            y_5 = y_7([1,3,4,6,7], :);
            return
        elseif (length(size(y_7)) == 3)
            y_5 = y_7([1,3,4,6,7], :, :);
        end
    end

    function y_3 = state_523(y_5)
        if (length(size(y_5)) < 3)
            y_3 = y_5([2,4,5], :);
            return
        elseif (length(size(y_5)) == 3)
            y_3 = y_5([2,4,5], :, :);
        end
    end

    function [ddy, ddshape] = second_der_3_tube(s, y)
        % 3 tubes; y = [p,R,theta2,theta3,u1z,u2z,u3z,mbx,mby]  
        % r = 1;
        [ka1_s, ka2_s, ka3_s] = ka_s(s);
        p = y(1:3); R = reshape(y(4:12),3,3);
        theta2 = y(13);  theta3 = y(14);
        u1z = y(15);  u2z = y(16);  u3z = y(17);
        mbx = y(18);  mby = y(19);
        % V = d/dxk[theta2,theta3,u1z,u2z,u3z,mbx,mby]
        mbz = G * (J1 * u1z + J2 * u2z + J3 * u3z);
        mb = [mbx; mby; mbz];
        sum_EIK = E * (I1*ka1_s*[1;0] + ...
                       I2*ka2_s*[cos(theta2);sin(theta2)] + ...
                       I3*ka3_s*[cos(theta3);sin(theta3)] );
        u1xy = ([mbx; mby] + sum_EIK) / (E * (I1 + I2 + I3));
        u1y = u1xy(2);
        u2y = [- sin(theta2), cos(theta2)] * u1xy;
        u3y = [- sin(theta3), cos(theta3)] * u1xy;
        u1 = [u1xy; u1z];

        e3 = [0;0;1];
        dp = R * e3;
        dR = R * hat(u1); dR_vec = reshape(dR, 9, 1);
        dtheta2 = u2z - u1z;
        dtheta3 = u3z - u1z;
        du1z = - (ka1_s * E * I1 / (G * J1)) * u1y;
        du2z = - (ka2_s * E * I2 / (G * J2)) * u2y;
        du3z = - (ka3_s * E * I3 / (G * J3)) * u3y;
        dmbxy = select_xy(- hat(u1) * mb - hat(e3)*R'*F);

        % V = d/dxk[theta2,theta3,u1z,u2z,u3z,mbx,mby]

        O3I3 = [zeros(3,3), eye(3)];
        VV = [dtheta2;dtheta3;du1z;du2z;du3z;dmbxy];
        EE = [dp; inv_hat(dR * R')];
        sum_dEIK = E * (...
           I2*ka2_s*[-sin(theta2);
                     cos(theta2)]*VV(1,:) + ...
           I3*ka3_s*[-sin(theta3)
                     cos(theta3)]*VV(2,:) );
        pu1xy_pxk = 1 / (E * (I1 + I2 + I3)) * (VV(6:7, :) + sum_dEIK);
        pu1_pxk = [pu1xy_pxk; VV(3, :)];
        pu1y_pxk = pu1xy_pxk(2, :);
        pu2y_pxk = VV(1,:) * ([-cos(theta2), - sin(theta2)] * u1xy) + ...
                    [-sin(theta2), cos(theta2)]*pu1xy_pxk;
        pu3y_pxk = VV(2,:) * ([-cos(theta3), - sin(theta3)] * u1xy) + ...
                    [-sin(theta3), cos(theta3)]*pu1xy_pxk;
        pmb_pxk = [VV(6:7, :); G*(J1*VV(3,:)+J2*VV(4,:)+J3*VV(5,:))];
        e3_dR1dxk_F = - hat(e3)*R'*cross([zeros(3),eye(3)]*EE,F);
        d_pmbxy_pxk_ds = [eye(2), zeros(2,1)] * ...
                     (- cross(pu1_pxk, mb) ...
                      - cross(u1, pmb_pxk) ...
                      - e3_dR1dxk_F);
        % creat dV=ddy_ds^2; dE=[ddp,ddR];
        dV = [
           VV(4,:) - VV(3,:);
           VV(5,:) - VV(3,:);
           - (ka1_s * E * I1) / (G * J1) * pu1y_pxk;
           - (ka2_s * E * I2) / (G * J2) * pu2y_pxk;
           - (ka3_s * E * I3) / (G * J3) * pu3y_pxk;
           d_pmbxy_pxk_ds];

        dE = [-hat(R*e3) * O3I3 * EE;
             R * pu1_pxk];
        
        ddy = dV;
        ddshape = dE;
    end

    function [ddy, ddshape] = second_der_2_tube(s, y)
        % 2 tubes; y = [p,R,theta2,u1z,u2z,mbx,mby]  
        % r = 1;
        [ka1_s, ka2_s, ka3_s] = ka_s(s);
        p = y(1:3); R = reshape(y(4:12),3,3);
        theta2 = y(13);  
        u1z = y(14);  u2z = y(15);  
        mbx = y(16);  mby = y(17);
        % V = d/dxk[theta2,u1z,u2z,mbx,mby]
        mbz = G * (J1 * u1z + J2 * u2z);
        mb = [mbx; mby; mbz];
        sum_EIK = E * (I1*ka1_s*[1;0] + ...
                       I2*ka2_s*[cos(theta2);sin(theta2)]);
        u1xy = ([mbx; mby] + sum_EIK) / (E * (I1 + I2));
        u1y = u1xy(2);
        u2y = [- sin(theta2), cos(theta2)] * u1xy;
        u1 = [u1xy; u1z];

        e3 = [0;0;1];
        dp = R * e3;
        dR = R * hat(u1); dR_vec = reshape(dR, 9, 1);
        dtheta2 = u2z - u1z;
        du1z = - (ka1_s * E * I1 / (G * J1)) * u1y;
        du2z = - (ka2_s * E * I2 / (G * J2)) * u2y;
        dmbxy = select_xy(- hat(u1) * mb - hat(e3)*R'*F);

        % V = d/dxk[theta2,u1z,u2z,mbx,mby]

        O3I3 = [zeros(3,3), eye(3)];
        VV = [dtheta2;du1z;du2z;dmbxy];
        EE = [dp; inv_hat(dR * R')];
        sum_dEIK = E * (I2*ka2_s*[-sin(theta2);
                     cos(theta2)]*VV(1,:));
        pu1xy_pxk = 1 / (E * (I1 + I2)) * (VV(4:5, :) + sum_dEIK);
        pu1_pxk = [pu1xy_pxk; VV(2, :)];
        pu1y_pxk = pu1xy_pxk(2, :);
        pu2y_pxk = VV(1,:) * ([-cos(theta2), - sin(theta2)] * u1xy) + ...
                    [-sin(theta2), cos(theta2)]*pu1xy_pxk;
        pmb_pxk = [VV(4:5, :); G*(J1*VV(2,:)+J2*VV(3,:))];
        e3_dR1dxk_F = - hat(e3)*R'*cross([zeros(3),eye(3)]*EE,F);
        d_pmbxy_pxk_ds = [eye(2), zeros(2,1)] * ...
                     (- cross(pu1_pxk, mb) ...
                      - cross(u1, pmb_pxk) ...
                      - e3_dR1dxk_F);
        % creat dV=ddy_ds^2; dE=[ddp,ddR];
        dV = [
           VV(3,:) - VV(2,:);
           - (ka1_s * E * I1) / (G * J1) * pu1y_pxk;
           - (ka2_s * E * I2) / (G * J2) * pu2y_pxk;
           d_pmbxy_pxk_ds];

        dE = [-hat(R*e3) * O3I3 * EE;
             R * pu1_pxk];
        
        ddy = dV;
        ddshape = dE;
    end

    function [ddy, ddshape] = second_der_1_tube(s, y)
        % 3 tubes; y = [p,R,u1z,mbx,mby]  
        [ka1_s, ka2_s, ka3_s] = ka_s(s);
        p = y(1:3); R = reshape(y(4:12),3,3);
        u1z = y(13);
        mbx = y(14);  mby = y(15);
        % V = d/dxk[theta2,theta3,u1z,u2z,u3z,mbx,mby]
        mbz = G * (J1 * u1z);
        mb = [mbx; mby; mbz];
        sum_EIK = E * (I1*ka1_s*[1;0]);
        u1xy = ([mbx; mby] + sum_EIK) / (E * I1);
        u1y = u1xy(2);
        u1 = [u1xy; u1z];

        e3 = [0;0;1];
        dp = R * e3;
        dR = R * hat(u1);
        du1z = - (ka1_s * E * I1 / (G * J1)) * u1y;
        dmbxy = select_xy(- hat(u1) * mb - hat(e3)*R'*F);

        % V = d/dxk[theta2,theta3,u1z,u2z,u3z,mbx,mby]

        O3I3 = [zeros(3,3), eye(3)];
        VV = [du1z;dmbxy];
        EE = [dp; inv_hat(dR * R')];
        sum_dEIK = 0;
        pu1xy_pxk = 1 / (E * I1) * (VV(2:3, :) + sum_dEIK);
        pu1_pxk = [pu1xy_pxk; VV(1, :)];
        pu1y_pxk = pu1xy_pxk(2, :);
        pmb_pxk = [VV(2:3, :); G*J1*VV(1,:)];
        e3_dR1dxk_F = - hat(e3)*R'*cross([zeros(3),eye(3)]*EE,F);
        d_pmbxy_pxk_ds = [eye(2), zeros(2,1)] * ...
                     (- cross(pu1_pxk, mb) ...
                      - cross(u1, pmb_pxk) ...
                      - e3_dR1dxk_F);
        % creat dV=ddy_ds^2; dE=[ddp,ddR];
        dV = [
           - (ka1_s * E * I1) / (G * J1) * pu1y_pxk;
           d_pmbxy_pxk_ds];

        dE = [-hat(R*e3) * O3I3 * EE;
             R * pu1_pxk];
        
        ddy = dV;
        ddshape = dE;
    end

    if nargin == 4
        s_in = ctr.L1 + q(2);
    end

    alpha1 = q(1); beta1 = q(2); alpha2 = q(3);
    beta2 = q(4); alpha3 = q(5); beta3 = q(6);
    F = w(1:3);  L = w(4:6);
    Ls3=ctr.Ls3; Ls2=ctr.Ls2; Ls1=ctr.Ls1;
    Lc3=ctr.Lc3; Lc2=ctr.Lc2; Lc1=ctr.Lc1;
    s_32 = Ls3 + Lc3 + beta3;  % transition 3 -> 2 tubes
    s_21 = Ls2 + Lc2 + beta2;  % transition 2 -> 1 tubes
    s_end = Ls1 + Lc1 + beta1;  % end potision
    E = ctr.E; G = ctr.G;
    ka1 = ctr.ka1; ka2 = ctr.ka2; ka3 = ctr.ka3;
    I1=ctr.I1; I2=ctr.I2; I3=ctr.I3; J1=ctr.J1; J2=ctr.J2; J3=ctr.J3;
    ode_solver = ctr.ode_solver;

    % yu = [u1z0, u2z0, u3z0, mbx0, mby0], get V0, E0 from yu
    yu_0 = yu0_true;
    u1z0 = yu_0(1); u2z0 = yu_0(2); u3z0 = yu_0(3);
    theta2_0 = alpha2 - alpha1 - (beta2 * u2z0 - beta1 * u1z0);
    theta3_0 = alpha3 - alpha1 - (beta3 * u3z0 - beta1 * u1z0);
    p0 = zeros(3,1);
    R1_0 = [cos(alpha1-beta1*u1z0), -sin(alpha1-beta1*u1z0), 0;
            sin(alpha1-beta1*u1z0), cos(alpha1-beta1*u1z0),  0;
            0,                      0,                       1];
    % note that E0 is HYBRID JACOBIAN form !!
    V0 = [-1,u1z0,1,-u2z0,0,    0,zeros(1,6),beta1,-beta2,     0,0,0;
          -1,u1z0,0,    0,1,-u3z0,zeros(1,6),beta1,     0,-beta3,0,0;
          zeros(5,12),                                        eye(5)];
    U0 = zeros(7,17,17); 
    U0(:,:,2) = [zeros(2,12), [1;1], zeros(2,4); zeros(5,17)];   %pV0_pbeta1
    U0(:,:,4) = [zeros(2,13), [-1;0], zeros(2,3); zeros(5,17)];  %pV0_pbeta2
    U0(:,:,6) = [zeros(2,14), [0;-1], zeros(2,2); zeros(5,17)];  % pV0_pbeta3
    U0(1:2,2,13) = [1;1];
    U0(1,4,14) = -1;
    U0(2,6,15) = -1;

    E0_c1 = [0,0,0,0,0,1]';
    E0_c2 = [0,0,0,0,0,-u1z0]';
    E0_c13 = [0,0,0,0,0,-beta1]';
    E0 = [E0_c1, E0_c2, zeros(6,10), E0_c13, zeros(6,4)];
    D0 = zeros(6,17,17);
    % D0(:,:,2) = [zeros(3,17); zeros(3,12), [0;0;-1], zeros(3,4)];  %pE0_pbeta1
    D0(6,13,2) = -1;
    D0(6,2,13) = -1;

    ini_3_tube = [p0; reshape(R1_0, 9,1); theta2_0; theta3_0;
        reshape(yu_0,5,1); reshape(V0,7*17,1); reshape(E0,6*17,1);
        reshape(U0,7*17*17,1); reshape(D0,6*17*17,1)];
    % disp('Calculate Hessian using 1 IVP: ')
    % Tube 3 -> Tube 2, gives y+|s=beta3;
    [s_s32, y_s32] = ode_solver(@ode_Hess_3_tube, [0, s_32], ini_3_tube);

    % ------------ 3 tubes shift to 2 tubes. -----------

    y_s3r_end = y_s32(end, :)';
    [p_s32,R_s32,state_7_s32,V_s3r,E_s3r,U_s3r,D_s3r] = unzip_3_tube(y_s3r_end);
    g_s3r = [R_s32, p_s32; zeros(1,3), 1];
    dy_s3r_end = ode_Hess_3_tube(s_32, y_s3r_end);
    [dp_s3r,dR_s3r,dy_7_s3r,dV_s3r,dE_s3r,~,~] = unzip_3_tube(dy_s3r_end);
    dg_s3r = [dR_s3r, dp_s3r; zeros(1,4)];
    [ddy_s3r, ddshape_s3r] = second_der_3_tube(s_32,y_s3r_end(1:19));

    y_s2left_pR = [p_s32;reshape(R_s32,9,1);state_725(state_7_s32)];
    dy_s2left_end = ode_Hess_2_tube(s_32, y_s2left_pR);
    [dp_s2left,dR_s2left,dy_5_s2left,~,~,~,~] = unzip_2_tube(dy_s2left_end);
    dg_s2left = [dR_s2left, dp_s2left; zeros(1,4)];
    [ddy_s2left, ddshape_s2left] = second_der_2_tube(s_s32, y_s2left_pR);
    % [ddy_s2left, ddshape_s2left] = second_der_2_tube(s_32,y_s2left_pR);

    V_s2left = state_725(V_s3r);
    V_s2left(:, 6) = V_s2left(:, 6) + state_725(dy_7_s3r) - dy_5_s2left;
    E_s2left = E_s3r;
    E_s2left(:,6) = E_s2left(:,6) + [dp_s3r;inv_hat(dR_s3r*R_s32')] ...
        - [dp_s2left;inv_hat(dR_s2left*R_s32')];
    
    y_s2left_VE = [y_s2left_pR; reshape(V_s2left,5*17,1); reshape(E_s2left,6*17,1)];
    dy_s2left_VE = ode_Hess_2_tube(s_32, y_s2left_VE);
    [~,~,~,dV_s2left,dE_s2left,~,~] = unzip_2_tube(dy_s2left_VE);

    % first order transition
    V_s32_full_der = V_s3r;
    V_s32_full_der(:,6) = V_s32_full_der(:,6) + dy_7_s3r;
    U_s32_full_der = U_s3r;
    U_s32_full_der(:,[1:5,7:17],6) = U_s32_full_der(:,[1:5,7:17],6) + ...
        dV_s3r(:,[1:5,7:17]);
    U_s32_full_der(:,6,[1:5,7:17]) = U_s32_full_der(:,6,[1:5,7:17]) + ...
        reshape(dV_s3r(:,[1:5,7:17]),[7,1,16]);
    
    U_s2left = state_725(U_s3r);
    U_s2left(:,[1:5,7:17],6) = U_s2left(:,[1:5,7:17], 6) + ...
        state_725(dV_s3r(:,[1:5,7:17])) - dV_s2left(:,[1:5,7:17]);
    U_s2left(:,6,[1:5,7:17]) = U_s2left(:,6,[1:5,7:17]) ...
        + reshape(state_725(dV_s3r(:,[1:5,7:17])),[5,1,16]) ...
        - reshape(dV_s2left(:,[1:5,7:17]),[5,1,16]);

    O3_I3 = [zeros(3,3), eye(3)];
    D_s2left = D_s3r;
    D_s2left(:, [1:5,7:17],6) = D_s2left(:, [1:5,7:17],6) + ...
        dE_s3r(:, [1:5,7:17]) - dE_s2left(:,[1:5,7:17]);
    for kk = [1:5,7:17]
        pR32_pxk = hat(O3_I3 * E_s3r(:,kk)) * R_s32;
        D_s2left(:,6,kk) = D_s2left(:,6,kk) + dE_s3r(:,kk) - dE_s2left(:,kk) ...
            +[zeros(3,1); inv_hat(-(dR_s2left-dR_s3r)*pR32_pxk'+pR32_pxk*(dR_s2left-dR_s3r)')];
    end

    D_s2left_t = D_s3r;
    D_s2left_t(:, [1:5,7:17],6) = D_s2left_t(:, [1:5,7:17],6) + ...
        dE_s3r(:, [1:5,7:17]);
    for kk = [1:5,7:17]
        pR32_pxk_t = hat(O3_I3 * E_s3r(:,kk)) * R_s32;
        D_s2left_t(:,6,kk) = D_s2left_t(:,6,kk) + dE_s3r(:,kk) ...
            +[zeros(3,1); inv_hat(-(-dR_s3r)*pR32_pxk_t'+pR32_pxk_t*(-dR_s3r)')];
    end


    % second order transition
    U_s32_full_der(:,6,6) = U_s32_full_der(:,6,6) + 2*dV_s3r(:,6) + ...
        ddy_s3r;

    U_s2left(:,6,6) = U_s2left(:,6,6) + 2*state_725(dV_s3r(:,6)) ...
        + state_725(ddy_s3r) - 2*dV_s2left(:,6) - ddy_s2left;

    pR_px6_s3r = hat(O3_I3 * E_s3r(:,6)) * R_s32;
    pR_px6_s2left = hat(O3_I3 * E_s2left(:,6)) * R_s32;
%     D_s2left(:,6,6) = D_s2left(:,6,6) + 2*dE_s3r(:,6) + ddshape_s3r ...
%         - 2*dE_s2left(:,6) - ddshape_s2left ...
%         + [zeros(3,1); inv_hat((pR_px6_s2left-pR_px6_s3r)*(pR_px6_s2left+dR_s2left)' ...
%                                 + (pR_px6_s2left+dR_s2left)*(dR_s2left-dR_s3r)')];
    D_s2left(:,6,6) = D_s2left(:,6,6) + 2*dE_s3r(:,6) + ddshape_s3r ...
        - 2*dE_s2left(:,6) - ddshape_s2left ...
        + [zeros(3,1);inv_hat(dR_s3r*pR_px6_s3r'-pR_px6_s3r*dR_s3r')] ...
        - [zeros(3,1);inv_hat(dR_s2left*pR_px6_s2left'-pR_px6_s2left*dR_s2left')];
    D_s2left_t(:,6,6) = D_s2left_t(:,6,6) + 2*dE_s3r(:,6) + ddshape_s3r ...
        + [zeros(3,1);inv_hat(dR_s3r*pR_px6_s3r'-pR_px6_s3r*dR_s3r')];
%     D_s2left(:,6,1:5) = D_s2left(:,6,1:5) + ...
%         reshape(dE_s3r(:,1:5),[6,1,5]) - ...
%         reshape(dE_s2left(:,1:5),[6,1,5]);
%     D_s2left(:,6,6);
    % U_s2left()

    % -------- 2 tube -> 1 tube ----------

    ini_2_tube = [y_s2left_VE; reshape(U_s2left,5*17*17,1);reshape(D_s2left,6*17*17,1)];
    [s_s21, y_s21] = ode_solver(@ode_Hess_2_tube, [s_32,s_21], ini_2_tube);
    y_s2r_end = y_s21(end,:)';
    [p_s21,R_s21,state_5_s21,V_s2r,E_s2r,U_s2r,D_s2r] = unzip_2_tube(y_s2r_end);
    dy_s2r_end = ode_Hess_2_tube(s_21, y_s2r_end);
    [dp_s2r,dR_s2r,dy_5_s2r,dV_s2r,dE_s2r,~,~] = unzip_2_tube(dy_s2r_end);
    [ddy_s2r, ddshape_s2r] = second_der_2_tube(s_21,y_s2r_end(1:17));

    y_s1left_pR = [p_s21;reshape(R_s21,9,1);state_523(state_5_s21)];
    dy_s1left_end = ode_Hess_1_tube(s_21, y_s1left_pR);
    [dp_s1left, dR_s1left, dy_3_s1left,~,~,~,~] = unzip_1_tube(dy_s1left_end);
    [ddy_s1left, ddshape_s1left] = second_der_1_tube(s_21,y_s1left_pR);

    V_s1left = state_523(V_s2r);
    V_s1left(:, 4) = V_s1left(:,4) + state_523(dy_5_s2r) - dy_3_s1left;
    E_s1left = E_s2r;
    E_s1left(:,4) = E_s1left(:,4) + [dp_s2r;inv_hat(dR_s2r*R_s21')] ...
                                - [dp_s1left;inv_hat(dR_s1left*R_s21')];

    y_s1left_VE = [y_s1left_pR;reshape(V_s1left,3*17,1);reshape(E_s1left,6*17,1)];
    dy_s1left_VE = ode_Hess_1_tube(s_21, y_s1left_VE);
    [~,~,~,dV_s1left, dE_s1left,~,~] = unzip_1_tube(dy_s1left_VE);

    % first order transition
    V_s21_full_der = V_s2r;
    V_s21_full_der(:,4) = V_s21_full_der(:,4) + dy_5_s2r;
    U_s21_full_der = U_s2r;
    U_s21_full_der(:,[1:3,5:17],4) = U_s21_full_der(:,[1:3,5:17],4) + ...
        dV_s2r(:, [1:3,5:17]);
    U_s21_full_der(:,4,[1:3,5:17]) = U_s21_full_der(:,4,[1:3,5:17]) + ...
        reshape(dV_s2r(:, [1:3,5:17]),[5,1,16]);

    U_s1left = state_523(U_s2r);
    U_s1left(:,[1:3,5:17],4) = U_s1left(:,[1:3,5:17],4) + ...
        state_523(dV_s2r(:,[1:3,5:17])) - dV_s1left(:,[1:3,5:17]);
    U_s1left(:,4,[1:3,5:17]) = U_s1left(:,4,[1:3,5:17]) ...
        + reshape(state_523(dV_s2r(:,[1:3,5:17])),[3,1,16]) ...
        - reshape(dV_s1left(:,[1:3,5:17]),[3,1,16]);

    D_s1left = D_s2r;
    D_s1left(:,[1:3,5:17],4) = D_s1left(:,[1:3,5:17],4) + dE_s2r(:,[1:3,5:17]) ...
        - dE_s1left(:,[1:3,5:17]);
    for kk = [1:3,5:17]
        pR21_pxk = hat(O3_I3 * E_s2r(:,kk)) * R_s21;
        D_s1left(:,4,kk) = D_s1left(:,4,kk) + dE_s2r(:,kk) - dE_s1left(:,kk) ...
            +[zeros(3,1); inv_hat(-(dR_s1left-dR_s2r)*pR21_pxk'+pR21_pxk*(dR_s1left-dR_s2r)')];
    end

    D_s1left_t = D_s2r;
    D_s1left_t(:,[1:3,5:17],4) = D_s1left_t(:,[1:3,5:17],4) + dE_s2r(:,[1:3,5:17]);
    for kk = [1:3,5:17]
        pR21_pxk = hat(O3_I3 * E_s2r(:,kk)) * R_s21;
        D_s1left_t(:,4,kk) = D_s1left_t(:,4,kk) + dE_s2r(:,kk)...
            +[zeros(3,1); inv_hat(-(0-dR_s2r)*pR21_pxk'+pR21_pxk*(0-dR_s2r)')];
    end

    % second oder transition
    U_s21_full_der(:,4,4) = U_s21_full_der(:,4,4) + 2*dV_s2r(:,4) + ddy_s2r;
    U_s1left(:,4,4) = U_s1left(:,4,4) + 2*state_523(dV_s2r(:,4)) ...
        + state_523(ddy_s2r) - 2*dV_s1left(:,4) - ddy_s1left;

    pR_px4_s2r = hat(O3_I3 * E_s2r(:,4)) * R_s21;
    pR_px4_s1left = hat(O3_I3 * E_s1left(:,4)) * R_s21;
    D_s1left(:,4,4) = D_s1left(:,4,4) + 2*dE_s2r(:,4) + ddshape_s2r ...
        -2*dE_s1left(:,4) - ddshape_s1left ...
        + [zeros(3,1);inv_hat(dR_s2r*pR_px4_s2r'-pR_px4_s2r*dR_s2r')] ...
        - [zeros(3,1);inv_hat(dR_s1left*pR_px4_s1left'-pR_px4_s1left*dR_s1left')];

    D_s1left_t(:,4,4) = D_s1left_t(:,4,4) + 2*dE_s2r(:,4) + ddshape_s2r ...
        + [zeros(3,1);inv_hat(dR_s2r*pR_px4_s2r'-pR_px4_s2r*dR_s2r')];

    % ------------ 1 tube -> end -----------
    ini_1_tube = [y_s1left_VE; reshape(U_s1left,3*17*17,1);reshape(D_s1left,6*17*17,1)];
    [s_s1e,y_s1e] = ode_solver(@ode_Hess_1_tube, [s_21, s_end], ini_1_tube);
    y_s1r_end = y_s1e(end, :)';
    [p_s1e,R_s1e,state_3_s1e,V_s1r,E_s1r,U_s1r,D_s1r] = unzip_1_tube(y_s1r_end);
    dy_s1r_end = ode_Hess_1_tube(s_end, y_s1r_end);
    [dp_s1r,dR_s1r,dstate_3_s1r,dV_s1r,dE_s1r,~,~] = unzip_1_tube(dy_s1r_end);
    [ddp_s1r, ddshape_s1r] = second_der_1_tube(s_end, y_s1r_end(1:15));    

    g_end = [R_s1e, p_s1e; zeros(1,3), 1];

    % transition of V, E
    V_end = V_s1r;
    V_end(:, 2) = V_end(:, 2) + dstate_3_s1r;
    E_end = E_s1r;
    E_end(:,2) = E_end(:,2) + [dp_s1r; inv_hat(dR_s1r*R_s1e')];

    % first order transition
    U_end = U_s1r;
    U_end(:,[1,3:17],2) = U_end(:,[1,3:17],2) + dV_s1r(:,[1,3:17]);
    U_end(:,2,[1,3:17]) = U_end(:,2,[1,3:17]) + reshape(dV_s1r(:,[1,3:17]),[3,1,16]);

    D_end = D_s1r;
    D_end(:,[1,3:17],2) = D_end(:,[1,3:17],2) + dE_s1r(:,[1,3:17]);
    for kk = [1,3:17]
        pR_1e_pxk = hat(O3_I3 * E_s1r(:,kk)) * R_s1e;
        D_end(:,2,kk) = D_end(:,2,kk) + dE_s1r(:,kk) ...
            +[zeros(3,1); inv_hat(-(-dR_s1r)*pR_1e_pxk'+pR_1e_pxk*(-dR_s1r)')];
    end

    % second order transition
    U_end(:,2,2) = U_end(:,2,2) + 2*dV_s1r(:,2) + ddp_s1r;

    pR_px2_s1r = hat(O3_I3 * E_s1r(:,2)) * R_s1e;
    D_end(:,2,2) = D_end(:,2,2) + 2*dE_s1r(:,2) + ddshape_s1r ...
        + [zeros(3,1);inv_hat(dR_s1r*pR_px2_s1r'-pR_px2_s1r*dR_s1r')];
    % find B, DB from U, V  only work when L=0 !!
    pL_pxk = [zeros(3,9),eye(3),zeros(3,5)];
    O3_I3 = [zeros(3,3), eye(3)];
    e3 = [0;0;1];
    p_RTL_pxk_end = R_s1e'*hat(L)*O3_I3*E_end + R_s1e'*pL_pxk;
    
    u1z_l1 = state_3_s1e(1);
    u2z_l2 = state_5_s21(3);
    u3z_l3 = state_7_s32(5);
    mbxy_l1 = state_3_s1e(2:3);

    b_end = [G * J1 * u1z_l1 - e3' * R_s1e' * L;
         G * J2 * u2z_l2; 
         G * J3 * u3z_l3;
         mbxy_l1 - select_xy(R_s1e' * L)];

    B_end = [
        G*J1*V_end(1,:) - e3'*p_RTL_pxk_end;
        G*J2*V_s21_full_der(3, :);
        G*J3*V_s32_full_der(5, :);
        V_end(2:3, :) - p_RTL_pxk_end(1:2, :);
    ];

    DB = zeros(5,17,17);
    pJ_px = zeros(6,6,17);
    pC_px = zeros(6,6,17);
    Bq = B_end(:, 1:6); Bw = B_end(:, 7:12); Bu = B_end(:,13:17); 
    Eq_end = E_end(:, 1:6); Ew_end = E_end(:, 7:12); Eu_end = E_end(:,13:17); 
    for rr = 1:17
        ppRTL_pxkpxr = R_s1e'*hat(L)*O3_I3*D_end(:,:,rr) ...
            + (-R_s1e'*hat(O3_I3*E_end(:,rr))*hat(L)+R_s1e'*hat(pL_pxk(:,rr)))*O3_I3*E_end ...
            - R_s1e'*hat(O3_I3*E_end(:,rr))*pL_pxk;
        DB(:,:,rr) = [G*J1*U_end(1,:,rr) - e3'*ppRTL_pxkpxr;
                      G*J2*U_s21_full_der(3,:,rr);
                      G*J3*U_s32_full_der(5,:,rr);
                      U_end(2:3,:,rr) - ppRTL_pxkpxr(1:2,:)];
        pJ_px(:,:,rr) = D_end(:,1:6,rr) - D_end(:,13:17,rr)*(Bu\Bq) ...
            + Eu_end*(Bu\DB(:,13:17,rr))*(Bu\Bq) ...
            - Eu_end*(Bu\DB(:,1:6,rr));
        pC_px(:,:,rr) = D_end(:,7:12,rr) - D_end(:,13:17,rr)*(Bu\Bw) ...
            + Eu_end*(Bu\DB(:,13:17,rr))*(Bu\Bw) ...
            - Eu_end*(Bu\DB(:,7:12,rr));
    end 

    % find the intermedia value

    if norm(s_in - s_end) < 1e-15
        E_in = E_end;
        D_in = D_end;
    elseif s_21 < s_in && s_in < s_end
        y_in = interp1(s_s1e, y_s1e, s_in);
        [p_in,R_in,~,V_in,E_in,U_in,D_in] = unzip_1_tube(y_in);
    elseif s_32 < s_in && s_in <= s_21
        y_in = interp1(s_s21, y_s21, s_in);
        [p_in,R_in,~,V_in,E_in,U_in,D_in] = unzip_2_tube(y_in);
    elseif 0 <= s_in && s_in <= s_32
        y_in = interp1(s_s32, y_s32, s_in);
        [p_in,R_in,~,V_in,E_in,U_in,D_in] = unzip_3_tube(y_in);
    else
        disp('no matched s_in')
    end

    Eq_in = E_in(:, 1:6); Ew_in = E_in(:, 7:12); Eu_in = E_in(:,13:17); 
    for rr = 1:17
        pJ_px_in(:,:,rr) = D_in(:,1:6,rr) - D_in(:,13:17,rr)*(Bu\Bq) ...
            + Eu_in*(Bu\DB(:,13:17,rr))*(Bu\Bq) ...
            - Eu_in*(Bu\DB(:,1:6,rr));
        pC_px_in(:,:,rr) = D_in(:,7:12,rr) - D_in(:,13:17,rr)*(Bu\Bw) ...
            + Eu_in*(Bu\DB(:,13:17,rr))*(Bu\Bw) ...
            - Eu_in*(Bu\DB(:,7:12,rr));
    end 

    % the first order staff
    J_end = Eq_end - Eu_end * inv(Bu) * Bq;
    C_end = Ew_end - Eu_end * inv(Bu) * Bw;

    J_in = Eq_in - Eu_in * inv(Bu) * Bq;
    C_in = Ew_in - Eu_in * inv(Bu) * Bw;
    % hessian is the jacobian of jacobian of shape.
    % we also need the jacobian of compliance of shape for redundancy
    % resolution
    Jq = pJ_px(:,:,1:6); Jw = pJ_px(:,:,7:12); Ju = pJ_px(:,:,13:end);
    Cq = pC_px(:,:,1:6); Cw = pC_px(:,:,7:12); Cu = pC_px(:,:,13:end);
    hess = Jq - tensorprod(Ju, inv(Bu)*Bq, 3, 1);
    jacob_compl = Cq - tensorprod(Cu, inv(Bu)*Bq, 3, 1);

    Jq_in = pJ_px_in(:,:,1:6); Jw_in = pJ_px_in(:,:,7:12); Ju_in = pJ_px_in(:,:,13:end);
    Cq_in = pC_px_in(:,:,1:6); Cw_in = pC_px_in(:,:,7:12); Cu_in = pC_px_in(:,:,13:end);
    hess_in = Jq_in - tensorprod(Ju_in, inv(Bu)*Bq, 3, 1);
    jacob_compl_in = Cq_in - tensorprod(Cu_in, inv(Bu)*Bq, 3, 1);

end
