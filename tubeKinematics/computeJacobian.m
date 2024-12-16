function J = computeJacobian(R,p)


    [~, ns] = size(p);

    p_skew = zeros(3,3,ns);

    p_skew(1,2,:) = -p(3,:);
    p_skew(2,1,:) = p(3,:);
    p_skew(1,3,:) = p(2,:);
    p_skew(3,1,:) = -p(2,:);
    p_skew(2,3,:) = -p(1,:);
    p_skew(3,2,:) = p(1,:);

    p_skew = reshape(p_skew, 3, []);

    % compute J
    % J = zeros(3*ns, 3*ns);
    % for i = 2:ns
    %     J(3*i-2:3*i, 1:3*i-3) = p_skew(:, 1:3*i-3) - repmat(p_skew(:, 3*i-2:3*i), 1, i-1);
    % end
    % 
    % for j = 1:ns
    %     J(:, 3*j-2:3*j) = J(:, 3*j-2:3*j) * R(:,:,j);
    % end

    J = repmat(p_skew, ns, 1);
    J = J + J';
    
    for j = 1:ns
        J(1:3*j, 3*j-2:3*j) = 0;
        J(:, 3*j-2:3*j) = J(:, 3*j-2:3*j) * R(:,:,j);
    end


end