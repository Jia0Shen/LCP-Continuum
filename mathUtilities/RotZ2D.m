function Rz = RotZ2D(theta)

st = sin(theta);
ct = cos(theta);
Rz = [ct -st; st ct];

end
