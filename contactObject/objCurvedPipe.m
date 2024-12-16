classdef objCurvedPipe
   properties
      g0   % 4x4xn 
      p0   % 3xn
      r
      mu
      T_history
      cornerFlag
   end
   methods
      function C = objCurvedPipe(p0, g0, r, mu)
            C.p0 = p0;
            C.r = r;
            C.mu = mu;
            C.g0 = g0;
            C.cornerFlag = false;
      end

      %% depth function
      function [ds, ns, p_proj] = computeDepth(C, p)

        p3n = reshape(p, 3, []);
        % np = size(p3n, 2);
        % np0 = size(C.p0, 2);
        
        idx = knnsearch(C.p0', p3n');
        p_proj_center = C.p0(:, idx);

        ds = C.r - vecnorm(p_proj_center-p3n, 2, 1);

        if nargout == 2
            ns = normalize(p_proj_center-p3n, 1, "norm");
        elseif nargout == 3
            ns = normalize(p_proj_center-p3n, 1, "norm");
            p_proj = p_proj_center - C.r * ns;
        end

    end

    %% projection function
    function p_proj = project(C, p)
        [~, ~, p_proj] = computeDepth(C, p);
    end

    %% normal vector function
    function ns = getNormal(C, p)
        [~, ns] = computeDepth(C, p);
    end

    %% geometry mesh function
    function [X, Y, Z] = getMesh(C)

        T0 = eye(4);
        res = size(C.g0,3);

        rin = C.r;

        for i = 1:res
            C.g0(:,:,i) = T0*C.g0(:,:,i);
        end

        % Compute bishop frames
        for i = 2:res
            w = LogSO3(C.g0(1:3,1:3,i-1)'*C.g0(1:3,1:3,i));
            C.g0(1:3,1:3,i) = C.g0(1:3,1:3,i)*LargeSO3([0,0,-w(3)]);
        end

        m = 50; % points on cross sectional circle

        X = zeros(res, m);
        Y = zeros(res, m);
        Z = zeros(res, m);

        t = linspace(0,2*pi,m);
        
        inCircle = [rin*cos(t); rin*sin(t); zeros(1,m);ones(1,m)];
        for i = 1:res

            inCross = C.g0(:,:,i)*inCircle;

            X(i,:) = inCross(1,:);
            Y(i,:) = inCross(2,:);
            Z(i,:) = inCross(3,:);
        end

        % Plot
        % surf(X,Y,Z, 'MeshStyle', 'both', 'LineStyle', 'none', 'FaceAlpha', alpha, 'FaceColor', color);

    end

    %% detect contact: linear combination of two corner points
    function obsContact = detectContact(C, tube, p, cornerFlag)

        % shrink the pipe
        C.r = C.r - tube.rout;

        % p should be 3xn: s=0-L
        p = reshape(p, 3, []);
        obsContact = [];

        [d, ns, p_proj] = C.computeDepth(p);
        IdxBodySurf = find(d < 0);
        dBodySurf = - d(IdxBodySurf);  % get penetration depth

        % remove the consective index
        [IdxBodySurfSelect, dBodySurfSelect] = removeConsective(IdxBodySurf, dBodySurf);

        for Idxp = IdxBodySurfSelect % IdxBodySurf

            contactSurfI.type = 'surfaceContact';
            contactSurfI.tube_point = p(:, Idxp);
            contactSurfI.tube_point_id = Idxp;
            contactSurfI.point = p_proj(:, Idxp);
            contactSurfI.normal = ns(:, Idxp);
            contactSurfI.penetrateDepth = - d(Idxp);  %signed dpeth

            obsContact = [obsContact, contactSurfI];
        end

        % change back the pipe 
        C.r = C.r + tube.rout;

    end

   end
end