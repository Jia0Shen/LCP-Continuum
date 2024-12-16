function T = LargeSE3(w,v)

wsqr = w'*w;
wnorm = sqrt(wsqr);

if wnorm > eps
    R = LargeSO3(w);
    W = [0 -w(3) w(2);w(3) 0 -w(1);-w(2) w(1) 0];
    Wv = W*v;

%     P = eye(3) + (1-cos(wnorm))/wsqr*W + (wnorm - sin(wnorm))/(wnorm*wsqr) * W^2;
%     T = [R P*v; zeros(1,3) 1];
    Pv = v + (1-cos(wnorm))/wsqr*Wv + (wnorm - sin(wnorm))/(wnorm*wsqr) * W * Wv;
    T = [R Pv; zeros(1,3) 1];
else
    T = [eye(3) v; zeros(1,3) 1];
end

% end

