function B = hat_mat(A)
    dim = size(A);  % A:3xn, or 3x1xn
    n = dim(end); % pick n from 3xn or 3x1xn;
    B = zeros(3,3,n);
    B(1,2,:) = - A(3,:); B(2,1,:) = A(3,:);
    B(1,3,:) = A(2,:);   B(3,1,:) = - A(2,:); 
    B(2,3,:) = - A(1,:); B(3,2,:) = A(1,:);
end