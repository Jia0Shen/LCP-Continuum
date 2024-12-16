function C = createPlaneEdge(p1, p2, n, mu)
% p1: first end point of the edge
% p2: second end point of the edge 
% z: normal direction of the plane
% mu: friction coefficient

C = struct;
C.p1 = p1;
C.p2 = p2;
C.n = n/norm(n);
C.mu = mu;

% functions
C.computeDepth = @(p)computeDepth(p);
C.project = @(p)project(p);

C.getMesh = @(num) getMesh(num);
C.detectContact = @(tube, p, cornerFlag)detectContact(tube, p, cornerFlag);

%% depth function
    function d = computeDepth(p)
        d = -n' * (p - p1);
    end

%% projection function
    function q = project(p)
        % q = p projected onto plane
        q = p - n*(n'*(p-p0));
    end

%% geometry mesh function
    function [X, Y, Z] = getPatch(width)
        if nargin == 0
            width = 3;   % 3mm width default
        end

        coord = [   0    0    0;
                    0.5  0    0;
                    0.5  0.5  0;
                    0    0.5  0;
                    0    0    0.5;
                    0.5  0    0.5;
                    0.5  0.5  0.5;
                    0    0.5  0.5;];
        idx = [4 8 5 1 4; 
            1 5 6 2 1; 
            2 6 7 3 2; 
            3 7 8 4 3; 
            5 8 7 6 5; 
            1 4 3 2 1]';
        xc = coord(:,1);
        yc = coord(:,2);
        zc = coord(:,3);

        X = xc(idx); Y = yc(idx); Z = zc(idx);

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
        
        % p_body = p(:, p(3,:) >= h);

        p_center = p0 + h/2*z;

        if isempty(p_body)
            IdxBody = [];
        else
            p_body_proj = (eye(3)-z*z')*(p_body-p_center);
            dist_body_proj = vecnorm(p_body_proj, 2, 1);
            dist_height = z' * (p_body - p_center);
            dist_body_corner = -tube.rout + sqrt((dist_body_proj-r).^2+dist_height.^2);
    
            % find the closest to the corner
            if cornerFlag
                cornerContactRange = 1;
            else
                cornerContactRange = 0;   % close to 1mm is corner contact
            end

            % IdxBody = find(D < cornerContactRange);
            [Dmin, IdxBody] = min(abs(dist_body_corner));
            % disp('Dmin is ' +string(Dmin))
    
            % dirTagent = normalize(p_body(:,IdxBody+1)-p_body(:,IdxBody-1), 'norm');
            dirTagent = normalize(p(:,IdxBody+1)-p(:,IdxBody), 'norm');
            dir_radius = normalize(p(:,IdxBody), 'norm');
            normProj = normalize((eye(3)-dirTagent*dirTagent')*dir_radius,'norm');

            if Dmin > cornerContactRange 
                IdxBody = [];
            end
    
            if z'*normProj < 0
                IdxBody = [];
            end

        end
        % IdxCornerContact = IdxCorner(IdxBody);

        IdxBodyRemove = [];
        
        if ~isempty(IdxBody)
            contactI.point = p_body_proj(:,IdxBody)+p_center;
            contactI.tube_point_id = IdxBody; 
            contactI.tube_point = p_body(:,IdxBody);
            contactI.type = 'cornerContact';
            contactI.penetrateDepth = dist_body_corner(IdxBody);

            % dirTagent = normalize(p_body(:,IdxBody+1)-p_body(:,IdxBody-1), 'norm');
            % normProj = - (eye(3)-dirTagent'*dirTagent)*(contactI.point-contactI.tube_point);
            % normProj = contactI.tube_point-contactI.point;
            % dir_radius = - normalize(contactI.point-p_center, 'norm');
            % dir_radius = normalize(contactI.point-p_center, 'norm');
            % normProj = normalize((eye(3)-dirTagent*dirTagent')*dir_radius,'norm');
            contactI.normal = normProj;

            obsContact = [obsContact, contactI];
        
            % 2. then find surface contacts
            removeRange = 5;  % remove the ovrelap within 5 discretization point
            IdxBodyRemove = IdxBody + (-removeRange:removeRange);
            IdxBodyRemove = unique(IdxBodyRemove(:));
            IdxBodyRemove(IdxBodyRemove>=size(p,2)) = [];  % never remove the tip.
        end
        % only care about the part of tube that inside the cylinder
        p_care = p(:, z'*(p-p_center)<=0);

        [d] = computeDepth(p_care, false);
        IdxBodySurf = find(d > - tube.rout);
        dBodySurf = - d(IdxBodySurf);

        % remove the consective index
        [IdxBodySurfSelect, dBodySurfSelect] = removeConsective(IdxBodySurf, dBodySurf);

        for Idxp = IdxBodySurfSelect % IdxBodySurf

            % check if it overlaps with corner
            if find(IdxBodyRemove==Idxp)
                continue
            else
                contactSurfI.type = 'surfaceContact';
                contactSurfI.tube_point = p_care(:, Idxp);
                contactSurfI.tube_point_id = Idxp;
                contactSurfI.point = project(p_care(:, Idxp));
                contactSurfI.normal = getNormal(p_care(:, Idxp));
                contactSurfI.penetrateDepth = d(Idxp);  %signed dpeth
            end

            obsContact = [obsContact, contactSurfI];
        end      

    end

end