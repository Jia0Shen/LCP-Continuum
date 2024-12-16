function C = createCapsule(p0, z, r, h, mu, in_or_out)
% p0: cylinder center
% z: axis of 
% r: radius
% h: height of central cylinder (excluding half spheres)
% mu: friction coefficient
% in_or_out: either 'in' or 'out' to specify the side of contact surface
% by default, in_or_out = 'out' to contact outer surface of cylinder


C = struct;
C.p0 = p0;
C.z = z/norm(z);
C.r = r;
C.h = h;
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

C.getMesh = @(temp)getMesh(temp);


%% depth function
    function d = computeDepth(p)
        dp = p - p0;

        pz = z'*dp;
        pz(pz > 0.5*h) = 0.5*h;
        pz(pz < -0.5*h) = -0.5*h;

        dp = dp - z*pz;
        d = r - sqrt(sum(dp.^2));

        d = C.normal_dir * d;
    end

%% projection function
    function q = project(p)
        dp = p - p0;

        pz = z'*dp;
        pz(pz > 0.5*h) = 0.5*h;
        pz(pz < -0.5*h) = -0.5*h;

        dp = dp - z*pz;

        dist = sqrt(sum(dp.^2));
        dist(dist < 1e-16) = 1e-16;
        q = p - dp +  dp ./ dist * r;
    end

%% normal vector function
    function n_ret = getNormal(p)
        dp = p - p0;

        pz = z'*dp;
        pz(pz > 0.5*h) = 0.5*h;
        pz(pz < -0.5*h) = -0.5*h;

        dp = dp - z*pz;

        dist = sqrt(sum(dp.^2));
        dist(dist < 1e-16) = 1e-16;
        n_ret = dp ./ dist;

        n_ret = C.normal_dir * n_ret;
    end

%% geometry mesh function
    function [X, Y, Z] = getMesh(temp)

        [a, b, c] = cylinder(1,50);
        
        c = (c - 0.5) * h/r;
        
        [aa, bb, cc] = sphere(50);
        
        id = [26:50, 1:26];
        aa = aa(:,id);
        bb = bb(:,id);
        cc = cc(:,id);
        
        a = [aa(1:25,:); a; aa(27:51,:)];
        b = [bb(1:25,:); b; bb(27:51,:)];
        c = [cc(1:25,:) - 0.5*h/r; c; cc(27:51,:) + 0.5*h/r];
        
        
        
        N = null(z');
        R = [N, z];
        
        XYZ = cell(1,3);
        for i = 1:3
            XYZ{i} = r * (R(i,1) * a + R(i,2) * b) + r * R(i,3) * c + p0(i);
        end
        X = XYZ{1};
        Y = XYZ{2};
        Z = XYZ{3};

    end

end