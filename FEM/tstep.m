function [a] = tstep(C, K, a0, f, edof, ex, ey)
%tstep function for time stepping
%   Solves implicit scheme to calculate evolution of transient temperature
%   distribution
a = 0;
t_tot = 600;
steps = 30;
dt = t_tot/steps;
aold = a0;
eT = extract(edof, aold);
    figure(101)
    hold on
    patch(ex',ey',eT', 'EdgeColor','none')
    patch(ex',-ey',eT', 'EdgeColor','none')
    colormap(hot);
    colorbar;  
    caxis([-60 35])
    xlabel('x-position [m]')
    ylabel('y-position [m]')
    title("Temperature distribution at t = 0")

for i = 1:steps
    anew = (C+dt*K)\(C*aold + dt*f);
    aold = anew;
    
    if i == round(steps/5)
        eT = extract(edof, aold);
        figure()
        hold on
        patch(ex',ey',eT', 'EdgeColor','none')
        patch(ex',-ey',eT', 'EdgeColor','none')
        colormap(hot);
        colorbar;  
        caxis([-60 35])
        xlabel('x-position [m]')
        ylabel('y-position [m]')
        title("Temperature distribution at t = " + i*dt)
        
    elseif i == round(steps/5*2)
        eT = extract(edof, aold);
        figure()
        hold on
        patch(ex',ey',eT', 'EdgeColor','none')
        patch(ex',-ey',eT', 'EdgeColor','none')
        colormap(hot);
        colorbar;  
        caxis([-60 35])
        xlabel('x-position [m]')
        ylabel('y-position [m]')
        title("Temperature distribution at t = " + i*dt)
        
    elseif i == round(steps/5*3)
        eT = extract(edof, aold);
        figure()
        hold on
        patch(ex',ey',eT', 'EdgeColor','none')
        patch(ex',-ey',eT', 'EdgeColor','none')
        colormap(hot);
        colorbar;  
        caxis([-60 35])
        xlabel('x-position [m]')
        ylabel('y-position [m]')
        title("Temperature distribution at t = " + i*dt)
        
    elseif i == round(steps/5*4)
        eT = extract(edof, aold);
        figure()
        hold on
        patch(ex',ey',eT', 'EdgeColor','none')
        patch(ex',-ey',eT', 'EdgeColor','none')
        colormap(hot);
        colorbar;
        caxis([-60 35])
        xlabel('x-position [m]')
        ylabel('y-position [m]')
        title("Temperature distribution at t = " + i*dt)
        
    end
end
    
    eT = extract(edof, aold);
    figure(102)
    hold on
    patch(ex',ey',eT', 'EdgeColor','none')
    patch(ex',-ey',eT', 'EdgeColor','none')
    colormap(hot);
    colorbar;
    caxis([-60 38])
    xlabel('x-position [m]')
    ylabel('y-position [m]')
    title("Temperature distribution at t = " + t_tot)

end

