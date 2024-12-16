function C = createCylinerCut(p0, z, r, h, mu, in_or_out)
% p0: cylinder center
% z: axis of 
% r: radius
% mu: friction coefficient
% h: height of central cylinder
% in_or_out: either 'in' or 'out' to specify the side of contact surface
% by default, in_or_out = 'out' to contact outer surface of cylinder

C = struct;
C.p0 = p0;
C.z = z/norm(z);
C.r = r;
C.mu = mu;
C.h = h;
C.in_or_out = in_or_out;

if strcmp(C.in_or_out, 'out')
    C.normal_dir = 1;
elseif strcmp(C.in_or_out, 'in')
    C.normal_dir = -1;
else
    error('Wrong input for in_or_out.');
end

% functions
C.computeDepth = @(p, flag)computeDepth(p, flag);
C.project = @(p)project(p);
C.getNormal = @(p)getNormal(p);

C.getMesh = @(num) getMesh(num);
C.detectContact = @(tube, p, cornerFlag)detectContact(tube, p, cornerFlag);

%% depth function
    function d = computeDepth(p, cutFlag)

    if nargin == 1
        cutFlag = true;
    end

        dp = p - p0;
        dp = dp - z*(z'*dp);
        d = r - sqrt(sum(dp.^2));

        d = C.normal_dir * d;

    if cutFlag
        p0z = p0(3);
        idx_cut = (p(3,:) > p0z+h/2 | p(3,:) < p0z-h/2);
        d(:, idx_cut) = nan;
    end

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
    function [X, Y, Z] = getMesh(edge_num)
        if nargin == 0
            edge_num = 50;
        end

        height = h;
        [a, b, c] = cylinder(1,edge_num);
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