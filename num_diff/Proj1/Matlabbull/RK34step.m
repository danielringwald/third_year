function [unew, err] = RK34step(f, told, uold, h)
unew = RK4step(f, told, uold, h); % 4th order approx

Y1=f(told,uold);
Y2=f(told + h/2, uold + h*Y1/2);
Y3=f(told + h/2, uold + h*Y2/2);
Z3=f(told + h, uold-h*Y1+2*h*Y2);
Y4=f(told + h,uold + h*Y3);

% Znew = uold + (h/6)*(Y1+4*Y2+Z3);  3rd order approx
err = norm((h/6)*(2*Y2 + Z3 - 2*Y3 -Y4)); 
end