function [shape1_t, R_trans, Errors] = compareShapes(shape1, shape2, s, rotFlag)
% input two shapes of 3xn. find min ||R*shape1-shape2|| if rotFlag is true.
    if nargin == 3
        % do not find rotation
        rotFlag = false;
    end

    len1 = size(shape1,2);  len2 = size(shape2,2);
    s1 = linspace(s(1),s(end),len1);
    s2 = linspace(s(1),s(end),len2);
    if len1 ~= length(s)
        shape1 = interp1(s1,shape1',s)';
    end

    if len2 ~= length(s)
        shape2 = interp1(s2,shape2',s)';
    end

    if ~rotFlag

        Errors.tip = norm(shape1(:,end) - shape2(:,end));
        Errors.meanBody = mean(vecnorm(shape1-shape2,2,1));
        Errors.maxBody = max(vecnorm(shape1-shape2,2,1));
        R_trans = [];
        shape1_t = shape1;
    else
        % use fminsearch/fmincon to find the rotation between two shape
        obj_fun = @(x) findRotErr(shape1, shape2, x);
        x0 = [0, 0, 0];
    
        [xsol, feval] = fminsearch(obj_fun, x0);
    
        R_trans = RotX(xsol(1))*RotY(xsol(2))*RotZ(xsol(3));
        shape1_t = R_trans *shape1;
        
        Errors.tip = norm(shape1_t(:,end) - shape2(:,end));
        Errors.meanBody = mean(vecnorm(shape1_t-shape2,2,1));
        Errors.maxBody = max(vecnorm(shape1_t-shape2,2,1));

    end

    function BodyErr = findRotErr(shape1, shape2, EulerAngles)
        R = RotX(EulerAngles(1))*RotY(EulerAngles(2))*RotZ(EulerAngles(3));
        % find min ||R*shape1-shape2||
        BodyErr = mean(vecnorm(R*shape1-shape2,2,1));
    end

        
end