function [tip_traj, whole_traj, contact_traj] = getPath2contact(tube, obstacles, traj, iniState)
    % inistate include the initial tube shape and env position.
    % find traj.
    alpha_traj = traj(1,:);
    beta_traj = traj(2,:);
    nt = size(traj,2);

    alpha0 = traj(1,1);

    % solve the shape from non-contact shape
    u = iniState.u;
    p = iniState.p;

    p0 = obstacles{2}.p0;
    obstacles{2}.T_history(:,:,1) = [eye(3), p0; 0 0 0 1];
    
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

        obstacles{2}.T_history(:,:,end+1) = [RotZ(alphaI), p0+[0;0;betaI]; 0 0 0 1];
        obstacles{2} = obstacles{2}.rebuild();

        [u, contacts, T, R, p] = getFrictionalContactShape3(tube, obstacles, u, contacts, CornerRange2);

        tip_traj(:,i) = p(:,end);
        whole_traj(:,i) = reshape(p,[],1);

    end

end