function [es] = therm2te(ny, E, alpha, dT)
%therm2te thermal stress for 2d triangle element
%  

es = alpha * E * dT / (1 - 2*ny) * [1 1 0]';

end

