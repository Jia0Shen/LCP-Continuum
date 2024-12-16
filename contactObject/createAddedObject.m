function C = createAddedObject(obs)
% obs: 1-d cell of objects
% mu is taken from the 1st object.

nObs = length(obs);

C = struct;
C.obs = obs;
C.nObs = nObs;

C.mu = obs{1}.mu;

% functions
C.computeDepth = @(p)computeDepth(p);
C.project = @(p)project(p);
C.getNormal = @(p)getNormal(p);

C.getMesh = @(temp)getMesh(temp);


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

end