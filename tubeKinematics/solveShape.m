function [T, R, p] = solveShape(T_base, u, s)


    n = length(s);

    T = zeros(4,4,n);

    T(:,:,1) = T_base;
    for i = 1:n-1

        ds = s(i+1) - s(i);
        
        % T(:,:,i+1) = T(:,:,i) * LargeSE3(u(:,i)*ds, [0;0;1]*ds);
        T(:,:,i+1) = T(:,:,i) * LargeSE3(u(:,i)*ds, [0;0;0]) * LargeSE3([0;0;0], [0;0;ds]);
    end


    % return R
    if nargout >= 2
        R = T(1:3,1:3,:);
    end


    % return p
    if nargout >= 3
        p = T(1:3,4,:);
        p = reshape(p,3,[]);
    end

end