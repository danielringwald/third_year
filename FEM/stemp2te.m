function [K_s, a, f_s] = stemp2te(K, edges_conv, p, ep, alpha_newt, dof, edof, Tinf, Tc)
%stemp2te Solver for stationary temperature conditions
%   K is global matrix
%   edges_conv is vctor with nodes on boundary
%   p are points of nodes
%   ep is thickness
%   alpha_newt: constant for Netown convection
%   dof: vector with all degrees of freedom
%   edof: connects elements with nodes
%   Tinf: left side temperature for convection
%   Tc: right side temperature for convetion

nen = 3;
K_s = K;
coord = p'/100;
f_s = zeros(size(p,2), 1);   %Preallocate memory to night f vector

for bn = 1:length(edges_conv(1,:))
    edge_nodes = edges_conv(:, bn);
    x1 = coord(edge_nodes(1), 1);
    y1 = coord(edge_nodes(1), 2);
    x2 = coord(edge_nodes(2), 1);
    y2 = coord(edge_nodes(2), 2);
    L = sqrt( (x2-x1)^2 + (y2-y1)^2);
    
    Kce = ep*alpha_newt*L*[1/3 1/6
                        1/6 1/3];
                    
    if x1 < max(coord(:,1))/2
        fce = ep*alpha_newt*Tinf*L*[1/2
                                 1/2];
    else
        fce = ep*alpha_newt*Tc*L*[1/2
                               1/2];
    end    
    
    indx = [edge_nodes(1) edge_nodes(2)];
    K_s(indx,indx) = K_s(indx,indx)+Kce;
    f_s(indx) = f_s(indx) + fce;
end

[ex, ey] = coordxtr(edof, coord, dof, nen);
a = solve(K_s, f_s);
eT = extract(edof, a);
figure()
hold on
patch(ex',ey',eT', 'EdgeColor','none')
patch(ex',-ey',eT', 'EdgeColor','none')

if Tinf == 40
    title('Stationary temperature distribution during day time[C]')
else
    title('Stationary temperature distribution during night time[C]')
end
colormap(hot);
colorbar;
xlabel('x-position [m]')
ylabel('y-position [m]')
axis equal
end

