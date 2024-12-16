function [shape1_t, R_trans, Errors] = compareMultiShapes(shape1, shape2, s, rotFlag)
% input two shapes of 3xn. find min ||R*shape1-shape2|| if rotFlag is true.
    if nargin == 3
        % do not find rotation
        rotFlag = false;
    end

    if ~iscell(shape1) || ~iscell(shape2)
        [shape1_t, R_trans, Errors] = compareShapes(shape1, shape2, s, rotFlag);
    else
        if ~rotFlag
            shape1_t = shape1;
            R_trans = eye(3);
            [~, Errors] = findRotMeanErr(shape1, shape2, [0,0,0], s);
        else
            % find the mean rotation
            obj_fun = @(x) findRotMeanErr(shape1, shape2, x, s);
            x0 = [0, 0, 0];

            % options = optimset('Display','iter', 'MaxIter', 100);
            % options = optimset('MaxIter', 200);
    
            % [xsol, feval] = fminsearch(obj_fun, x0, options);
            [xsol, feval] = fminsearch(obj_fun, x0);

            R_trans = RotX(xsol(1))*RotY(xsol(2))*RotZ(xsol(3));

            % [~, Errors] = findRotMeanErr(shape1, shape2, xsol, s);

            shape1_t = cell(size(shape1));

            for jj =  1:length(shape1)
                shape1_t{jj} = R_trans * shape1{jj};
            end

            % find the error after rotation
            [~, Errors] = findRotMeanErr(shape1_t, shape2, x0, s);

        end

    end


    function [errMetric, Error_t] = findRotMeanErr(shape1, shape2, EulerAngles, s)
        R = RotX(EulerAngles(1))*RotY(EulerAngles(2))*RotZ(EulerAngles(3));
        tipErr= [];  
        bodyErr= [];  
        maxErr= [];  
        for j = 1:length(shape1)
            shape1J = shape1{j};  shape2J = shape2{j};

            [~, ~, ErrorsJ] = compareShapes(R*shape1J, shape2J, s, false);
            tipErr = [tipErr, ErrorsJ.tip];
            bodyErr = [bodyErr, ErrorsJ.meanBody];
            maxErr = [maxErr, ErrorsJ.maxBody];
        end

        Error_t.tip = mean(tipErr);
        Error_t.meanBody = mean(bodyErr);
        Error_t.maxBody = mean(maxErr);

        errMetric = 1/3*(2*mean(bodyErr)+mean(tipErr));   % use meanbody as metric.
        % errMetric = mean(bodyErr);
    end
        
end