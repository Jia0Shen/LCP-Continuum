function plane = createPlane(p0, n, mu, alpha,dalpha)
% p0: point on plane
% n: normal vector
% mu: friction coefficient

plane = struct;
plane.p0 = p0;
plane.n = n / norm(n);
plane.mu = mu;
plane.T_history = [eye(3), p0; 0 0 0 1];
plane.cornerFlag = false;

% functions
plane.computeDepth = @(p)computeDepth(p);
plane.project = @(p)project(p);
plane.getNormal = @(p)getNormal(p);

plane.getMesh = @(width)getMesh(width);
plane.detectContact = @(tube, p, cornerFlag)detectContact(tube, p);

%% depth function
    function d = computeDepth(p)
        d = -n' * (p - p0);
    end

%% projection function
    function q = project(p)
        % q = p projected onto plane
        q = p - n*(n'*(p-p0));
    end

%% normal vector function
    function n_ret = getNormal(p)
        np = size(p,2);
        n_ret = n * ones(1,np);
    end

%% geometry mesh function
    function [X, Y, Z] = getMesh(width, shape)

        if (nargin < 1)
            width = 180;
        elseif (nargin == 1)
            shape = 'disk';          
        end
    
        % plane basis
        B = null(n');
    
        % plane corner points
        [X_mesh, Y_mesh] = meshgrid([-1,1] * width * 0.5, [-1,1] * width * 0.5);

        % if strcmp(shape, 'disk')
        %     circleMask = (X_mesh).^2 + (Y_mesh).^2 <= width^2/4;
        %     X_mesh(~circleMask) = NaN;
        %     Y_mesh(~circleMask) = NaN;
        % end
            
        p = reshape([X_mesh(:), Y_mesh(:)]*B', 2,2,3);

        X = p(:,:,1) + p0(1);
        Y = p(:,:,2) + p0(2);
        Z = p(:,:,3) + p0(3);
    end

%% detect contact
    function obsContact = detectContact(tube, p)

        % p should be 3xn: s=0-L
        obsContact = [];

        % Find surface contacts

        d = computeDepth(p);
        IdxBodySurf = find(d > - tube.rout);
        dBodySurf = - d(IdxBodySurf);

        % remove the consective index
        [IdxBodySurfSelect, dBodySurfSelect] = removeConsective(IdxBodySurf, dBodySurf);
        
        contactSurfI = [];
        for Idxp = IdxBodySurfSelect % IdxBodySurf
            contactSurfI.type = 'surfaceContact';
            contactSurfI.tube_point = p(:, Idxp);
            contactSurfI.tube_point_id = Idxp;
            contactSurfI.point = project(p(:, Idxp));
            contactSurfI.normal = getNormal(p(:, Idxp));
            contactSurfI.penetrateDepth = d(Idxp)+tube.rout;  %signed dpeth
        end

        obsContact = [obsContact, contactSurfI];
    end      

end