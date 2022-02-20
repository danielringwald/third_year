%Task 2.1
% An explicit Runge-Kutta method for IVP y'=f(t,y)
% of order 4, also called RK4
function unew = RK4step(f, told, uold, h)
Y1=f(told,uold); 
Y2=f(told + h/2, uold + h*Y1/2);
Y3=f(told + h/2, uold + h*Y2/2);
Y4=f(told + h,uold + h*Y3);

unew = uold + (h/6)*(Y1+2*Y2+2*Y3+Y4);
end