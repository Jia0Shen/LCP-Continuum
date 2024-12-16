classdef objPlane
   properties
      p0
      n
      mu
      T_history
      cornerFlag
   end
   methods
       function plane = objPlane(p0, n, mu)
           plane.p0 = p0;
           plane.n = n;
           plane.mu = mu;
           plane.cornerFlag = false;
       end

       function d = computeDepth(plane, p)
           d = - plane.n' * (p - plane.p0);
       end

       function plane = rebuild(plane)
           plane.p0 = plane.T_history(1:3,4,end);
       end

       function q = project(plane, p)
           % q = p projected onto plane
           q = p - plane.n*(plane.n'*(p-plane.p0));
       end

       function n_ret = getNormal(plane, p)
           np = size(p,2);
           n_ret = plane.n * ones(1,np);
       end

       function [X, Y, Z2] = getMesh(plane, radius, grids)
            if (nargin < 1)
                radius = 60;
            elseif nargin == 1
                grids = 100;
            end

            % plane basis
            % B = null(plane.n');

            width = 5;  % 5mm width plane.

            [X,Y,Z] = cylinder(radius,grids);
            Z2 = Z*width + plane.p0(3);

        end

        function obsContact = detectContact(plane, tube, p, cornerRange)

            % p should be 3xn: s=0-L
            obsContact = [];

            % Find surface contacts

            d = plane.computeDepth(p);
            IdxBodySurf = find(d > - tube.rout);
            dBodySurf = - d(IdxBodySurf);

            % remove the consective index
            [IdxBodySurfSelect, dBodySurfSelect] = removeConsective(IdxBodySurf, dBodySurf);

            contactSurfI = [];
            for Idxp = IdxBodySurfSelect % IdxBodySurf
                contactSurfI.type = 'surfaceContact';
                contactSurfI.tube_point = p(:, Idxp);
                contactSurfI.tube_point_id = Idxp;
                contactSurfI.point = plane.project(p(:, Idxp));
                contactSurfI.normal = plane.getNormal(p(:, Idxp));
                contactSurfI.penetrateDepth = d(Idxp)+tube.rout;  %signed dpeth
            end

            obsContact = [obsContact, contactSurfI];
        end


   end
end