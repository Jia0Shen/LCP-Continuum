function tube = CreatTube(ExposedLength, s_in, u_hat)
    
    if nargin == 1
        Mi = 400;
        len_precurve = 213.88;  % measurement
        s_pre = linspace(0, len_precurve, Mi);
        s_in = s_pre(s_pre > (len_precurve - ExposedLength)) - (len_precurve - ExposedLength);
        s_in(1) = 0;
        fitted_precurv = load("initial_uarray_cut_new.mat");
        u_hat_raw = fitted_precurv.u_array;
        uy_hat_raw = nonzeros(u_hat_raw);
        uy_hat = uy_hat_raw(s_pre > (len_precurve - ExposedLength));
        u_hat = uy2u(uy_hat);
    end
    
    tube = struct;

    tube.rin = 0.91/2;
    tube.rout = 1.32/2;  

    tube.v = 0.3;  
    tube.kb = 20.07e4;   % 20.07e-2 Nm2;
    
    tube.kt = tube.kb/(1+tube.v);
    
    tube.T_base = eye(4);
    
    % precurvature

    tube.s = s_in;
    tube.uhat = reshape(u_hat, 3, []);

end