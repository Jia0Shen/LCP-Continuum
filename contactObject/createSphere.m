function S = createSphere(p0, r, mu)
% p0: sphere center
% r: radius
% mu: friction coefficient


S = struct;
S.p0 = p0;
S.r = r;
S.mu = mu;

% functions
S.computeDepth = @(p)computeDepth(p);
S.project = @(p)project(p);
S.getNormal = @(p)getNormal(p);

S.getMesh = @getMesh;


%% depth function
    function d = computeDepth(p)
        dp = p - p0;
        d = r - sqrt(sum(dp.^2));
    end

%% projection function
    function q = project(p)
        dp = p - p0;
        dist = sqrt(sum(dp.^2));
        dist(dist < 1e-16) = 1e-16;
        q = (dp ./ dist) * r + p0;
    end

%% normal vector function
    function n_ret = getNormal(p)
        dp = p - p0;
        dist = sqrt(sum(dp.^2));
        dist(dist < 1e-16) = 1e-16;
        n_ret = dp ./ dist;
    end

%% geometry mesh function
    function [X, Y, Z] = getMesh()

        [X, Y, Z] = sphere();
        
        X = r*X + p0(1);
        Y = r*Y + p0(2);
        Z = r*Z + p0(3);

    end

end