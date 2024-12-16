function tube = changeTubeK(tube, bendScale, v)

    v0 = 0.3;

    if nargin == 1
        v = v0;
        bendScale = 1;
    elseif nargin == 2
        v = v0;
    end
    
    % scale = kb_straight / kb_curved; normally > 1
    % tube.v = v;
    % tube.kt = tube.kb / (1 + tube.v);

    tube.v = v0;

    n = length(tube.s);
    % curved arc length is 75 mm for our tube
    % len_curve = 75;  % mm
    % s_st = tube.s(end) - len_curve;
    s_st = 26;  % rigid before s = 26mm

    scaled_one = ones(1,n);
    scaled_one(tube.s < s_st) = bendScale;

    scaled_poission = 1 / (1 + v0) * ones(1,n);
    scaled_poission(tube.s <= s_st) = 1 / (1 + v);

    Kxy = [tube.kb; tube.kb] * scaled_one;

    Kz = tube.kb * scaled_one .* scaled_poission;

    tube.K = [Kxy; Kz];

end