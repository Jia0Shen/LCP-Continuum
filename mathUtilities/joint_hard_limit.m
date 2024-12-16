function q_clamped = joint_hard_limit(ctr, q)

        % joint hard limit 
        theta = (ctr.M_qt * q + ctr.N_qt);
        % gamma1 = theta(2); gamma2 = theta(4); gamma3 = theta(6);
        if (theta(2) < ctr.gamma1_lim(1) || theta(2) > ctr.gamma1_lim(2))
            disp('gamma1 hit limits.')
            theta(2) = max(theta(2),ctr.gamma1_lim(1));
            theta(2) = min(theta(2),ctr.gamma1_lim(2));
            % break;
        end

        if (theta(4) < ctr.gamma2_lim(1) || theta(4) > ctr.gamma1_lim(2))
            disp('gamma2 hit limits')
            theta(4) = max(theta(4),ctr.gamma2_lim(1));
            theta(4) = min(theta(4),ctr.gamma2_lim(2));
            % break;
        end

        if (theta(6) < ctr.gamma3_lim(1) || theta(6) > ctr.gamma3_lim(2))
            disp('gamma3 hit limits')
            theta(6) = max(theta(6),ctr.gamma3_lim(1));
            theta(6) = min(theta(6),ctr.gamma3_lim(2));
            % break;
        end

        q_clamped = (ctr.M_qt \ (theta - ctr.N_qt));
end