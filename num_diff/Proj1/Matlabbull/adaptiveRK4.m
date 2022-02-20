function [t,y] = adaptiveRK4(f, t0, tf, y0, tol)
t = [t0];
y = [y0];
k = 4;
h0ld = abs(tf-t0)*tol^(1/4) / (100*(1+norm(f(0,y0))));
errold = tol;

while t(end) < tf
    
    [unew, err] = RK34step(f, t(end), y(end), h0ld);
    hnew= newstep(tol, err, errold, h0ld, k);
    errold = err;
    h0ld = hnew;
    t(end + 1) = t(end) + hnew;
    y(end + 1) = unew; 
end

% subtrahera felet
[unew, err]     = RK34step(f, t(end-1), y(end-1), tf-t(end-1));
    y(end)          = unew;
    t(end)          = tf;
end