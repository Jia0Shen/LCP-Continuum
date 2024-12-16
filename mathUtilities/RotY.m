function Ry = RotY(theta)

st = sin(theta);
ct = cos(theta);
Ry = [ct  0  st;
      0   1   0;
      -st 0   ct];

end

