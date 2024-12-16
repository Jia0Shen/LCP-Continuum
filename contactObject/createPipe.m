function C = createPipe(obs)
% obs: 1-d cell of objects
% mu is taken from the 1st object.

nObs = length(obs);

C = struct;
C.obs = obs;
C.nObs = nObs;

C.mu = obs{1}.mu;
C.cornerFlag = true;
C.T_history = eye(4);

% functions
C.computeDepth = @(p)computeDepth(p);
C.project = @(p)project(p);
C.getNormal = @(p)getNormal(p);
C.getMesh = @(temp)getMesh(temp);
C.detectContact = @(tube, p, flag)detectContact(tube, p, flag);

%% depth function
    function [d, ids_obs] = computeDepth(p)
        depth = zeros(nObs, size(p,2));
        for io = 1:nObs
            depth(io,:) = obs{io}.computeDepth(p);            
        end

        [d, ids_obs] = min(depth);
    end

%% projection function
    function q = project(p)

        q = zeros(size(p));

        [~, ids_obs] = computeDepth(p);

        for io = 1:nObs
            cur_id = ids_obs == io;

            q(:,cur_id) = obs{io}.project(p(:,cur_id));       
        end
    end

%% normal vector function
    function n_ret = getNormal(p)

        n_ret = zeros(size(p));

        [~, ids_obs] = computeDepth(p);

        for io = 1:nObs
            cur_id = ids_obs == io;

            n_ret(:,cur_id) = obs{io}.getNormal(p(:,cur_id));       
        end
    end

%% geometry mesh function
    function [X, Y, Z] = getMesh(temp)
        for io = 1:nObs
            [X,Y,Z] = obs{io}.getMesh(temp);
            surf(X,Y,Z, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
            % surf(X,Y,Z, 'FaceAlpha', 0.5);
        end
    end

%% detect contact
    function obsContact = detectContact(tube, p, cornerFlag)
        if nargin == 2
            cornerFlag = true;  %if we model contact at corners
        end
        % p should be 3xn: s=0-L
        obsContact = [];

        % 1. first detect corner contact
        p_body = p(:, 1:end-1);  % tip cannot have corner contact
        
        z_mid = obs{1}.h;
        y_corner = 0; % linspace(-obs{1}.r, obs{1}.r, 100);
        x_corner = sqrt(obs{1}.r^2-y_corner.^2);
        z_corner = z_mid - sqrt(obs{2}.r^2 - y_corner.^2);
        corners = [x_corner; y_corner; z_corner];

        % find the closest to the corner
        % tic   
        [IdxCorner, D] = knnsearch(corners', p_body');  % D has columns as p
        % toc
        if cornerFlag
            cornerContactRange = 2;
        else
            cornerContactRange = 0;   % close to 2mm is corner contact
        end
        % IdxBody = find(D < cornerContactRange);
        [Dmin, IdxBody] = min(D);
        % disp('Dmin is ' +string(Dmin))
        if Dmin > cornerContactRange
            IdxBody = [];
        end
        % IdxCornerContact = IdxCorner(IdxBody);
        
        if ~isempty(IdxBody)
            for IdxBodyI = reshape(IdxBody, 1, [])
                IdxCornerContactI = IdxCorner(IdxBodyI);
                contactI.point = corners(:, IdxCornerContactI);
                contactI.tube_point_id = IdxBodyI; 
                contactI.tube_point = p_body(:,IdxBodyI);
                contactI.type = 'cornerContact';
                contactI.penetrateDepth = Dmin;
    
                dirTagent = normalize(p_body(:,IdxBodyI+1)-p_body(:,IdxBodyI-1), 'norm');
                % normProj = - (eye(3)-dirTagent'*dirTagent)*(contactI.point-contactI.tube_point);
                % normProj = contactI.tube_point-contactI.point;
                normProj = cross(dirTagent, [0,1,0]');
                % determine the normal direct
                % if contactI.tube_point(1) > obs{1}.r
                %     contactI.normal = normalize(normProj, 'norm');
                % else
                %     contactI.normal = normalize(normProj, 'norm');
                % end
                contactI.normal = normProj;
                

                obsContact = [obsContact, contactI];
            end
        end
        % 2. then find surface contacts
        removeRange = 25;  % remove the ovrelap within 5 discretization point
        IdxBodyRemove = IdxBody(:) + (-removeRange:removeRange);
        IdxBodyRemove = unique(IdxBodyRemove(:));
        IdxBodyRemove(IdxBodyRemove>=size(p,2)) = [];  % never remove the tip.

        [d, ids_obs] = computeDepth(p);
        IdxBodySurf = find(d > - tube.rout);
        dBodySurf = - d(IdxBodySurf);

        % remove the consective index
        [IdxBodySurfSelect, dBodySurfSelect] = removeConsective(IdxBodySurf, dBodySurf);

        for Idxp = IdxBodySurfSelect % IdxBodySurf

            % check if it overlaps with corner
            if find(IdxBodyRemove==Idxp)
                continue
            else
                idx_obsI = ids_obs(Idxp);
                contactSurfI.type = 'surfaceContact';
                contactSurfI.tube_point = p(:, Idxp);
                contactSurfI.tube_point_id = Idxp;
                contactSurfI.point = obs{idx_obsI}.project(p(:, Idxp));
                contactSurfI.normal = obs{idx_obsI}.getNormal(p(:, Idxp));
                contactSurfI.penetrateDepth = d(Idxp);  %signed dpeth
            end

            obsContact = [obsContact, contactSurfI];
        end      

    end
end