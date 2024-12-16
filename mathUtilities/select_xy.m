function A = select_xy(B)
    % B = [v1,v2,..,vn]; A = [v1xy, v2xy, ..., vnxy];
    A = B(1:2, :, :);
end