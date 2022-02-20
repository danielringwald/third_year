function [D] = isoDte(ny, E)
%isoDte creates constitutive matrix for isotropic material

D = E / ( (1+ny)*(1-2*ny) ) * [1-ny ny 0
                              ny 1-ny 0
                              0 0 1/2*(1-2*ny)];

end

