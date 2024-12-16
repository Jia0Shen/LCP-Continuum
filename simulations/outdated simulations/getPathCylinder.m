function [tip_traj, whole_traj, contact_traj] = getPathCylinder(tube, obstacles, traj, iniState, param)
    % inistate include the initial tube shape and env position.
    % find traj.
    if nargin == 4
        param = [];
    end
    
    alpha_traj = traj(1,:);
    beta_traj = traj(2,:);
    nt = size(traj,2);

    alpha0 = traj(1,1);

    % solve the shape from non-contact shape
    u = iniState.u;
    p = iniState.p;

    h1 = obstacles{1}.h;  h2 = obstacles{2}.h;
    h1_ini = 120;
    p1_ini = [0; 0; h1_ini-h1/2];
    p2_ini = [0; 0; h1_ini-h1-h2/2];
    p1 = obstacles{1}.p0;
    p2 = obstacles{2}.p0;
    p1p2 = {p1_ini, p2_ini};

    obstacles{1}.T_history(:,:,1) = [RotZ(alpha0), p1; 0 0 0 1];
    obstacles{1} = obstacles{1}.rebuild();

    obstacles{2}.T_history(:,:,1) = [RotZ(alpha0), p2; 0 0 0 1];
    obstacles{2} = obstacles{2}.rebuild();
    
    CornerRange1 = 0.5;
    CornerRange2 = 0.5;

    % ini
    tip_traj = zeros(3,nt);
    whole_traj = zeros(3*size(p,2), nt);
    contact_traj = cell(1,nt);

    tip_traj(:,1) = p(:,end);
    whole_traj(:,1) = reshape(p,[],1);
    
    for i = 2:nt
    
        contacts = [];
        [contacts, ~, ~, ~] = detectAdditionalContacts(tube, u, obstacles, contacts, CornerRange1);
        contact_traj{i-1} = contacts;

        alphaI = alpha_traj(i);
        betaI = beta_traj(i);

        obstacles{1}.T_history(:,:,end+1) = [RotZ(alphaI), p1+[0;0;betaI]; 0 0 0 1];
        obstacles{1} = obstacles{1}.rebuild();
    
        obstacles{2}.T_history(:,:,end+1) = [RotZ(alphaI), p2+[0;0;betaI]; 0 0 0 1];
        obstacles{2} = obstacles{2}.rebuild();
    
        % obstacles{1}.mu = 0.01;
        if ~isempty(param)
            [u, contacts, T, R, p] = getFrictionalContactShape3(tube, obstacles, u, contacts, CornerRange2,param);
            
        else
            [u, contacts, T, R, p] = getFrictionalContactShape3(tube, obstacles, u, contacts, CornerRange2);
        end
        tip_traj(:,i) = p(:,end);
        whole_traj(:,i) = reshape(p,[],1);

    end

end