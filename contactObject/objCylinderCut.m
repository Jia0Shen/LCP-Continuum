classdef objCylinderCut
   properties
      p0
      z
      r
      mu
      h
      in_or_out
      normal_dir
      T_history
      cornerFlag
   end
   methods
      function C = objCylinderCut(p0, z, r, h, mu, in_or_out)
            C.p0 = p0;
            C.z = z/norm(z);
            C.r = r;
            C.mu = mu;
            C.h = h;
            C.in_or_out = in_or_out;
            C.cornerFlag = true;  % default
            
            if strcmp(C.in_or_out, 'out')
                C.normal_dir = 1;
            elseif strcmp(C.in_or_out, 'in')
                C.normal_dir = -1;
            else
                error('Wrong input for in_or_out.');
            end

      end

      %% rebuilt
    function C = rebuild(C)
        C.p0 = C.T_history(1:3,4,end);
    end

      %% depth function
    function d = computeDepth(C, p, cutFlag)

        if nargin == 1
            cutFlag = true;
        end
    
            dp = p - C.p0;
            dp = dp - C.z*(C.z'*dp);
            d = C.r - sqrt(sum(dp.^2));
    
            d = C.normal_dir * d;
    
        if cutFlag
            p0z = C.p0(3);
            idx_cut = (p(3,:) > p0z+C.h/2 | p(3,:) < p0z-C.h/2);
            d(:, idx_cut) = nan;
        end

    end

          %% depth function
    function d = computeCutHeight(C, p)
        d = - p(3,:) + C.h/2 + C.p0(3);
    end

    %% projection function
    function q = project(C, p)
        dp = p - C.p0;
        dp = dp - C.z*(C.z'*dp);
        dist = sqrt(sum(dp.^2));
        dist(dist < 1e-16) = 1e-16;
        q = p - dp +  dp ./ dist * C.r;
    end

    %% project to corner function
    function q = project2corner(C,p,threshold)
        if nargin == 2
            threshold = 0;
        end

        q_surf = C.project(p);
        p_top = C.p0 + C.h/2*C.z;

        d_proj = dot(p_top-q_surf, C.z);

        q = q_surf + d_proj*C.z;

        % if d_proj < threshold
        %     % no project to corner
        %     q = q_surf;
        % else
        %     q = q_surf + (d_proj-threshold)*C.z;
        % end
    end

    %% normal vector function
    function n_ret = getNormal(C, p)

        dp = p - C.p0;
        dp = dp - C.z*(C.z'*dp);
        dist = sqrt(sum(dp.^2));
        dist(dist < 1e-16) = 1e-16;
        n_ret = dp ./ dist;
    
        n_ret = C.normal_dir * n_ret;

    end

%% geometry mesh function
    function [X, Y, Z] = getMesh(C, edge_num)
        if nargin == 0
            edge_num = 50;
        end

        [a, b, c] = cylinder(1,edge_num);
        c = c - 0.5;
        
        N = null(C.z');
        R = [N, C.z];

        XYZ = cell(1,3);
        for i = 1:3
            XYZ{i} = C.r * (R(i,1) * a + R(i,2) * b) + C.h * R(i,3) * c + C.p0(i);
        end
        X = XYZ{1};
        Y = XYZ{2};
        Z = XYZ{3};

    end

%% detect contact: single point
function obsContact = detectContactSingle(C, tube, p, cornerContactRange)
        if nargin == 2
            cornerContactRange = 0;  %if we model contact at corners
        end
        % p should be 3xn: s=0-L
        obsContact = [];

        % 1. first detect corner contact
        p_body = p(:, 1:end-1);  % tip cannot have corner contact
        
        % p_body = p(:, p(3,:) >= h);

        p_center = C.p0 + C.h/2*C.z;
        p_buttom = C.p0 - C.h/2*C.z;

        if isempty(p_body)
            IdxBody = [];
        else
            p_body_proj = (eye(3)-C.z*C.z')*(p_body-p_center);
            dist_body_proj = vecnorm(p_body_proj, 2, 1);
            dist_height = C.z' * (p_body - p_center);
            dist_body_corner = -tube.rout + sqrt((dist_body_proj-C.r).^2+dist_height.^2);
    
            % find the closest to the corner
            % if cornerFlag
            %     cornerContactRange = 0.5;
            % else
            %     cornerContactRange = 0;   % close to 0.5mm is corner contact
            % end

            % IdxBody = find(D < cornerContactRange);
            [Dmin, IdxBody] = min(abs(dist_body_corner));
            % disp('Dmin is ' +string(Dmin))
    
            % dirTagent = normalize(p_body(:,IdxBody+1)-p_body(:,IdxBody-1), 'norm');
            dirTagent = normalize(p(:,IdxBody+1)-p(:,IdxBody), 'norm');
            % dir_radius = normalize(p(:,IdxBody), 'norm');
            body_point = C.r*normalize(p_body_proj(:,IdxBody),'norm')+p_center;
            dir_radius = normalize(p_center-body_point, 'norm');
            % normProj = normalize((eye(3)-dirTagent*dirTagent')*dir_radius,'norm');
            rad_tang = dot(dirTagent, dir_radius);  % <0
            rad_mag = 1/sqrt(1+rad_tang^2);
            tang_mag = - rad_tang/sqrt(1+rad_tang^2);              
            normProj = rad_mag*dir_radius + tang_mag*dirTagent;

            if Dmin > cornerContactRange 
                IdxBody = [];
            end
    
            if C.z'*normProj < 0
                IdxBody = [];
            end

        end
        % IdxCornerContact = IdxCorner(IdxBody);

        IdxBodyRemove = [];
        
        if ~isempty(IdxBody)
            contactI.point = body_point;
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
        % p_care = p(:, C.z'*(p-p_center)<=0  &  C.z'*(p-p_buttom)>=0);
        logic_care = C.z'*(p-p_center)<=0  &  C.z'*(p-p_buttom)>=0;
        idx_p = 1:size(p,2);
        idx_care = idx_p(:, logic_care);
        p_care = p(:, logic_care);
        % p_care = p(:, C.z'*(p-p_center)<=0);

        [d] = C.computeDepth(p_care, false);
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
                contactSurfI.tube_point_id = idx_care(Idxp);
                contactSurfI.point = C.project(p_care(:, Idxp));
                contactSurfI.normal = C.getNormal(p_care(:, Idxp));
                contactSurfI.penetrateDepth = d(Idxp);  %signed dpeth
            end

            obsContact = [obsContact, contactSurfI];
        end      

end

%% detect contact: linear combination of two corner points
function obsContact = detectContact(C, tube, p, cornerContactRange)
        if nargin == 2
            cornerContactRange = 0;  %if we model contact at corners
        end
        % p should be 3xn: s=0-L
        obsContact = [];

        % 1. first detect corner contact
        p_body = p(:, 1:end-1);  % tip cannot have corner contact
        p_center = C.p0 + C.h/2*C.z;
        p_buttom = C.p0 - C.h/2*C.z;

        p_body_proj = (eye(3)-C.z*C.z')*(p_body-p_center);
        dist_body_proj = vecnorm(p_body_proj, 2, 1);
        dist_height = C.z' * (p_body - p_center);
        % dist_body_corner = -tube.rout + sqrt((dist_body_proj-C.r).^2+dist_height.^2);
        dist_body_corner = sqrt((dist_body_proj-C.r).^2+dist_height.^2);

        % find the closest to the corner
        d_pene = C.computeDepth(p_body, false);
        idxPene = find(d_pene(1:end-1)<0 & d_pene(2:end)>0);
        IdxCorner = [];

        for ip = idxPene
            ppi_1 = p_body(:,ip);  ppi = p_body(:,ip+1);
            % find the point on the cylinder from ppi_1 -> ppi:
            di_1 = abs(d_pene(ip));  di = abs(d_pene(ip+1));
            ppc = di/(di_1+di)*ppi_1 + di_1/(di_1+di)*ppi;
            dirTagentI_1 = normalize(ppi-ppi_1, 'norm');
            L = tube.rout / sqrt(1-(dirTagentI_1'*C.z)^2);
            % ppcRadiusOffest = ppc - L*C.z;
            ppcRadiusOffest = ppc;  
            
            % if ppc_z 
            if ppcRadiusOffest(3) < C.p0(3) + C.h/2 + cornerContactRange & ...
                  ppcRadiusOffest(3) > C.p0(3) + C.h/2 - cornerContactRange
                % corner contact detected
                dir_radius = [normalize(-ppcRadiusOffest(1:2), 'norm'); 0];
                rad_tang = dot(dirTagentI_1, dir_radius);  % <0
                % rad_mag = 1;
                % tang_mag = - rad_tang;              
                normProj = normalize(dir_radius - rad_tang*dirTagentI_1, 'norm');

                % only 1 corner contact?
                IdxCorner = ip;
                break;
            end
        end                

        % [Dmin, IdxBody] = mink(dist_body_corner, 2, 'ComparisonMethod', 'abs');

        IdxBodyRemove = [];
        
        if ~isempty(IdxCorner)
            % point is on surface (accurat), tube_point is the record of
            % previous discretized point.
            contactI.point = ppcRadiusOffest;
            contactI.tube_point_id = [IdxCorner,IdxCorner+1]; 
            contactI.tube_point = p_body(:,[IdxCorner,IdxCorner+1]);
            contactI.type = 'cornerContact';
            contactI.penetrateDepth = C.p0(3) + C.h/2 - ppcRadiusOffest(3);
            contactI.normal = normProj;

            obsContact = [obsContact, contactI];
        
            % 2. then find surface contacts
            removeRange = 5;  % remove the ovrelap within 5 discretization point
            IdxBodyRemove = IdxCorner + (-removeRange:removeRange);
            IdxBodyRemove = unique(IdxBodyRemove(:));
            IdxBodyRemove(IdxBodyRemove>=size(p,2)) = [];  % never remove the tip.
        end

        % only care about the part of tube that inside the cylinder
        logic_care = C.z'*(p-p_center)<=0  &  C.z'*(p-p_buttom)>=0;
        idx_p = 1:size(p,2);
        idx_care = idx_p(:, logic_care);
        p_care = p(:, logic_care);
        % p_care = p(:, C.z'*(p-p_center)<=0);

        [d] = C.computeDepth(p_care, false);
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
                contactSurfI.tube_point_id = idx_care(Idxp);
                contactSurfI.point = C.project(p_care(:, Idxp));
                contactSurfI.normal = C.getNormal(p_care(:, Idxp));
                contactSurfI.penetrateDepth = d(Idxp);  %signed dpeth
            end

            obsContact = [obsContact, contactSurfI];
        end      

    end

   end
end