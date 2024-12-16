function [contacts, T, R, p] = detectAdditionalContacts(tube, u, obstacles, contacts, cornerRange)

    % solve for shape
    [T, R, p] = solveShape(tube.T_base, u, tube.s);
    

    % compute depth w.r.t. obstacles
    nObs = length(obstacles);
    % ns = length(tube.s);

    % depth = zeros(nObs, ns);
    % for io = 1:nObs
    %     depth(io, :) = obstacles{io}.detectContact(tube, p);
    % end
    % 
    % % get points of penetration
    % ip = depth > 0; % points of penetration

    % compare with current contacts and add if necessary
    id_tol = 10;  % the contact within 10 discretization points are the same.
    isIDclose = @(id1, id2) abs(id1-id2) <= id_tol;
    for io = 1:nObs
        if obstacles{io}.cornerFlag
            cornerRangeI = cornerRange;
        else 
            cornerRangeI = 0;
        end
        obs_contact = obstacles{io}.detectContact(tube, p, cornerRangeI);
        % ip_cur = find(ip(io,:));
        for is = 1:length(obs_contact)
            % check existence
            bExist = false;
            for ic = 1:length(contacts)
                % if (contacts(ic).obstacle_id == io  ...
                %         && isIDclose(contacts(ic).tube_point_id,obs_contact(is).tube_point_id))
                if (isIDclose(contacts(ic).tube_point_id,obs_contact(is).tube_point_id))
                    bExist = true;
                end
            end
    
            % add if not exists
            if ~bExist
                % contacts(end+1) = contacts(end);
                % contacts(end).force(:) = 0;
                % contacts(end).obstacle_id = io;
                % contacts(end).tube_point_id = is;
                % contacts(end).point = obstacles{io}.project(p(:,is));
    
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%
                c = struct;
    
                c.force(:) = 0;
                c.obstacle_id = io;
                c.tube_point_id = obs_contact(is).tube_point_id;
                c.point = obs_contact(is).point;
                c.normal = obs_contact(is).normal;
                c.type = obs_contact(is).type;
                c.tube_point = obs_contact(is).tube_point;
                c.penetrateDepth = obs_contact(is).penetrateDepth;
    
                if isempty(contacts)
                    contacts = c;
                else
                    contacts(end+1) = c;
                end
            end
        end


    end


end