function C = createCyliner(p0, z, r, mu, in_or_out)
% p0: cylinder center
% z: axis of 
% r: radius
% mu: friction coefficient
% in_or_out: either 'in' or 'out' to specify the side of contact surface
% by default, in_or_out = 'out' to contact outer surface of cylinder


C = struct;
C.p0 = p0;
C.z = z/norm(z);
C.r = r;
C.mu = mu;
C.in_or_out = in_or_out;
if strcmp(C.in_or_out, 'out')
    C.normal_dir = 1;
elseif strcmp(C.in_or_out, 'in')
    C.normal_dir = -1;
else
    error('Wrong input for in_or_out.');
end

% functions
C.computeDepth = @(p)computeDepth(p);
C.project = @(p)project(p);
C.getNormal = @(p)getNormal(p);

C.getMesh = @(height)getMesh(height);
C.detectContact = @(tube, p, cornerFlag)detectContact(tube, p, cornerFlag);

%% depth function
    function d = computeDepth(p)
        dp = p - p0;
        dp = dp - z*(z'*dp);
        d = r - sqrt(sum(dp.^2));

        d = C.normal_dir * d;
    end

%% projection function
    function q = project(p)
        dp = p - p0;
        dp = dp - z*(z'*dp);
        dist = sqrt(sum(dp.^2));
        dist(dist < 1e-16) = 1e-16;
        q = p - dp +  dp ./ dist * r;
    end

%% normal vector function
    function n_ret = getNormal(p)
        dp = p - p0;
        dp = dp - z*(z'*dp);
        dist = sqrt(sum(dp.^2));
        dist(dist < 1e-16) = 1e-16;
        n_ret = dp ./ dist;

        n_ret = C.normal_dir * n_ret;
    end

%% geometry mesh function
    function [X, Y, Z] = getMesh(height)

        if nargin < 1
            height = 100;
        end

        [a, b, c] = cylinder(1,50);
        c = c - 0.5;
        
        N = null(z');
        R = [N, z];

        XYZ = cell(1,3);
        for i = 1:3
            XYZ{i} = r * (R(i,1) * a + R(i,2) * b) + height * R(i,3) * c + p0(i);
        end
        X = XYZ{1};
        Y = XYZ{2};
        Z = XYZ{3};

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