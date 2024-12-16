function ctr = build_ctr()
    % define the parameters needed for the CTR.
    ctr.E = 60e9;     % Pa. E = 2G(1+v)
    ctr.G = 23.1e9;   % Pa
    ctr.ka1 = 4; ctr.ka2 = 4; ctr.ka3 = 2;
    ctr.Ls1 = 450e-3; ctr.Lc1 = 150e-3; 
    ctr.L1 = ctr.Ls1 + ctr.Lc1;
    ctr.Ls2 = 250e-3; ctr.Lc2 = 150e-3; ctr.L2 = ctr.Ls2 + ctr.Lc2;
    ctr.Ls3 = 100e-3; ctr.Lc3 = 100e-3; ctr.L3 = ctr.Ls3 + ctr.Lc3;
    D1 = 1e-3; d1 = 0.5e-3;
    D2 = 1.75e-3; d2 = 1.25e-3;
    D3 = 2.5e-3; d3 = 2e-3;
    ctr.I1 = 1/64 * pi  * (D1^4 - d1^4);
    ctr.J1 = 2 * ctr.I1;
    ctr.I2 = 1/64 * pi  * (D2^4 - d2^4);
    ctr.J2 = 2 * ctr.I2;
    ctr.I3 = 1/64 * pi  * (D3^4 - d3^4);
    ctr.J3 = 2 * ctr.I3;
    ctr.K1 = diag([ctr.E*ctr.I1, ctr.E*ctr.I1, ctr.G*ctr.J1]);
    ctr.K2 = diag([ctr.E*ctr.I2, ctr.E*ctr.I2, ctr.G*ctr.J2]);
    ctr.K3 = diag([ctr.E*ctr.I3, ctr.E*ctr.I3, ctr.G*ctr.J3]);
    % transformation from q to theta: theta' = M*q' + N;
    ctr.M_qt =  [1 0 0 0 0 0; 
             0 1 0 -1 0 0; 
             0 0 1 0 0 0; 
             0 0 0 1 0 -1;
             0 0 0 0 1 0;
             0 0 0 0 0 1];
    ctr.N_qt = [0; ctr.L1-ctr.L2; 0; ctr.L2-ctr.L3; 0; ctr.L3];
    % joint limits
    ctr.gamma1_lim = [0.01, ctr.L1-ctr.L2];
    ctr.gamma2_lim = [0.01, ctr.L2-ctr.L3];
    ctr.gamma3_lim = [0.01, 0.2];

    % build the solver
    ctr.ode_solver = @ode45;
end