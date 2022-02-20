tol = 1e-8;
t0 = 0;
tf = 100;
%f = lotka();
y0 = [1; 1];

%a = 3;
%b = 9;
%dudt = @(t, u) [a*u(1), b*u(2)];
[t, y] = adaptiveRK4(lotka(), t0, tf, y0, tol)
