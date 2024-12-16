function adg = Adj(g)
    % g = [R, p; 0,0,0,1];
    R = g(1:3, 1:3); p = g(1:3, 4);
    adg = [R, hat(p) * R; 
           zeros(3), R];    
end