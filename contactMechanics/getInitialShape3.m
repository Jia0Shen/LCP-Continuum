function [u, contacts, T, R, p] = getInitialShape3(tube, obstacles, cornerFlag)
% Solve shape for tube-obstacle contacts
% constacts: array of struct(contact point, plane index, tube contact point index)
% obstacles: cells of obstacles

if nargin == 2
    cornerFlag = true;
end

    bPlot = false;

    nObs = length(obstacles);
    ns = length(tube.s);
    
    K = getTubeK(tube);
    invK = 1./K;
    
    % init variables
    f = zeros(nObs, ns);       % plane contact forces
    u = tube.uhat;


    % Check initial contact free configuration
    % [contacts, T, R, p] = detectAdditionalContacts(tube, u, obstacles, []);
    % if isempty(contacts)
    %     return;
    % end

    % loop for solving tube shape
    for iter = 1:100
        % calculate tube shape
        [~, R, p] = solveShape(tube.T_base, u, tube.s);
        if (bPlot && mod(iter-1,1) == 0)
            % Plot
            if (iter == 1)
                figure
                h = plot3Mat(p);
                hold on;
                for io = 1:nObs
                    plot_obj = obstacles{io};
                    [X,Y,Z] = plot_obj.getMesh(200);
                    % surf(X,Y,Z, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
                end

                axis equal
                grid on

                xlabel('x')
                ylabel('y')
            else
                set(h, 'XData', p(1,:), 'YData', p(2,:), 'ZData', p(3,:))
            end
            pause(0.2)
            drawnow

        end

        % find penetrations
        depth = zeros(nObs, ns);
        for io = 1:nObs
            depth(io, :) = obstacles{io}.computeDepth(p);
        end

%         if prod(prod(depth <= 1e-6))
%             break;
%         end

        % update f (force magnitudes)
        f = f + depth * 2e-8 * tube.kb;
        f(f < 0) = 0;

        % update u (J du ds = dp, J^T f = m = k(u-uhat))
        J = computeJacobian(R,p);

        f_sum = zeros(3,ns);    % force vectors
        for io = 1:nObs
            f_sum = f_sum + obstacles{io}.getNormal(p) .* f(io, :);
        end
 
        Jtf = J' * f_sum(:);
        u = reshape(invK .* Jtf, 3, []) + tube.uhat;

    end

    % % loop for solving tube shape
    % for iter = 1:100
    %     % calculate tube shape
    %     [~, R, p] = solveShape(tube.T_base, u, tube.s);
    %     if (bPlot && mod(iter-1,1) == 0)
    %         % Plot
    %         if (iter == 1)
    %             figure
    %             h = plot3Mat(p);
    %             hold on;
    %             for io = 1:nObs
    %                 [X,Y,Z] = obstacles{io}.getMesh(200);
    %                 surf(X,Y,Z);
    %             end
    % 
    %             axis equal
    %             grid on
    % 
    %             xlabel('x')
    %             ylabel('y')
    %         else
    %             set(h, 'XData', p(1,:), 'YData', p(2,:), 'ZData', p(3,:))
    %         end
    %         drawnow
    % 
    %     end
    % 
    %     % find penetrations
    %     depth = zeros(nObs, ns);
    %     for io = 1:nObs
    %         depth(io, :) = obstacles{io}.computeDepth(p);
    %     end
    % 
    %     if prod(prod(depth <= 1e-3))
    %         break;
    %     end
    % 
    % 
    %     % find f to yield depth = 0 for point currently with depth < 0
    %     J = computeJacobian(R,p);
    % 
    %     invK_m = u - tube.uhat;
    %     invK_m = invK_m(:);
    %     ds = tube.s(end) - tube.s(end-1); % assumed constant
    % 
    %     ip = depth > 0;
    % 
    %     for io = 1:nObs
    %         ip_cur = ip(io,:);
    % 
    %         % ip3{io} = repmat(ip_cur, 1, 3);
    %         % ip3{io} = ip3{io}(:)';
    % 
    %         normal{io} = obstacles{io}.getNormal(p(:, ip_cur));
    %         normal_all_s = zeros(3, ns);
    %         normal_all_s(:, ip_cur) = normal{io};
    % 
    %         N{io} = zeros(3*ns, ns);
    %         for k = 1:3
    %             N{io}(k:3*ns+3:end) = normal_all_s(k,:);
    %         end
    %         N{io} = N{io}(:, ip_cur);
    % 
    %         % Jp{io} = J(id{io}, id{io});
    %         % invKp{io} = invK(id{io})
    %         % 
    %         % invKp_mp{io} = invK_m(:, pen_id(io,:));
    %         % invKp_mp{io} = invKp_mp{io}(:);
    %         % 
    %         % A{io} = (Jp{io} .* invKp{io}) * Jp{io}' * ds;
    %         % q{io} = -Jp{io} * invKp_mp{io} * ds;
    %         % 
    %         % N{io} = 
    %         % 
    %         % nAn{io}
    % 
    % 
    %     end
    % 
    %     ip_unified = ip(1,:);
    %     N_unified = N{1};
    %     for io = 2:nObs
    %         ip_unified = ip_unified | ip(io,:);
    %         N_unified = [N_unified, N{io}];
    %     end
    % 
    %     ip3_unified = repmat(ip_unified, 3, 1);
    %     ip3_unified = ip3_unified(:);
    % 
    %     N_unified = N_unified(ip3_unified, :);
    % 
    %     Jp = J(ip3_unified, :);
    %     % invKp = invK(ip3_unified);
    % 
    %     % invKp_mp = invK_m(ip3_unified);
    % 
    %     % A = (Jp .* invKp') * Jp' * ds;
    %     % q = -Jp * invKp_mp * ds;
    % 
    %     A = (Jp .* invK') * Jp' * ds;
    %     q = -Jp * invK_m * ds;
    % 
    % 
    %     NtAN = N_unified'*A*N_unified;
    %     Ntq = N_unified'*q;
    % 
    %     % depth = NtAN * f + Ntq
    %     dd = depth(:);
    %     depth_unified = dd(ip(:));
    %     % f_unified = pinv(NtAN) * (depth_unified - Ntq);
    % 
    %     [~, f_unified] = LCPSolve(NtAN, Ntq - depth_unified);
    % 
    % 
    %     % update f (force magnitudes)
    % 
    %     sumf = abs(sum(f_unified));
    %     sumf_pos = sum(f_unified(f_unified > 0));
    % 
    %     f_unified(f_unified < 0) = 0;
    %     f_unified = f_unified * (sumf / sumf_pos);
    % 
    %     f(:) = 0;
    %     f(ip) = f_unified;
    %     f(f < 0) = 0;
    % 
    %     % update u (J du ds = dp, J^T f = m = k(u-uhat))
    %     J = computeJacobian(R,p);
    % 
    %     f_sum = zeros(3,ns);    % force vectors
    %     for io = 1:nObs
    %         f_sum = f_sum + obstacles{io}.getNormal(p) .* f(io, :);
    %     end
    % 
    %     Jtf = J' * f_sum(:);
    %     u = reshape(invK .* Jtf, 3, []) + tube.uhat;
    % 
    % 
    % end
    

    % constacts: array of struct(contact point, obstacle index, tube contact point index, force)
    contacts = [];
    for io = 1:nObs
        cid = find(f(io,:) > 0);

        for j = 1:length(cid)
            contact = struct;

            contact.point = obstacles{io}.project( p(:,cid(j)) );
            contact.obstacle_id = io;
            contact.force = obstacles{io}.getNormal(contact.point) * f(io, cid(j));
            contact.tube_point_id = cid(j);

            % add additional properties.
            contact.normal = obstacles{io}.getNormal( p(:,cid(j)) );
            contact.type = 'Undefined';
            contact.tube_point = p(:,cid(j));
            contact.penetrateDepth = - obstacles{io}.computeDepth( p(:,cid(j)) );

            % add to contacts
            if isempty(contacts)
                contacts = contact;
            else
                contacts(end+1) = contact;
            end
        end
    end

    for iter = 1:7
        contacts = [];
        [contacts, T, R, p] = detectAdditionalContacts(tube, u, obstacles, contacts, cornerFlag);
        [u, contacts, T, R, p] = getContactShape3(tube, obstacles, u, contacts, cornerFlag);
    end

    % T, R, p
    if nargout > 2
        [T, R, p] = solveShape(tube.T_base, u, tube.s);
    end



end