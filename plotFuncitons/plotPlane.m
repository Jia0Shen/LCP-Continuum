function h = plotPlane(plane, width)

    if (nargin < 2)
        width = 100;
    end
    
    % plane basis
    B = null(plane.n');

    % plane corner points
    [X,Y] = meshgrid([-1,1] * width * 0.5, [-1,1] * width * 0.5);
        
    p = reshape([X(:), Y(:)]*B', 2,2,3);
    h = surf(p(:,:,1) + plane.p0(1), p(:,:,2) + plane.p0(2), p(:,:,3) + plane.p0(3));

end