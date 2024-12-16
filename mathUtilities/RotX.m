function Rx = RotX(theta)

st = sin(theta);
ct = cos(theta);
% Rz = [ct -st 0; st ct 0; 0 0 1];
Rx = [1  0  0;
      0  ct  -st;
      0  st  ct];

end

