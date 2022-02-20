function dudt = lotka()
a = 3;
b = 9;
c = 15;
d = 15;
%u = [1, 1];
dudt = @(t, u) [a*u(1), b*u(2)];
end