function Rz = RotZ(theta)

st = sin(theta);
ct = cos(theta);
Rz = [ct -st 0; st ct 0; 0 0 1];

end

