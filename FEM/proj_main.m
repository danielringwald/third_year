%% Setup constans and geometry
%e är edges 
%p is points (coordinates)
%t triangles ([nodes], subdomain)

E_Ti = 110e9 ;     %Young's modulus [N / m^2]
E_Gl = 67e9 ;
ny_Ti = 0.34;           %Posisson's ratio [-]
ny_Gl = 0.2;
alpha_Ti = 9.4 * 10^-6;     %Expansions coefficient [1/K]
alpha_Gl = 7 * 10^-6;
rho_Ti = 4620 ;          %density [kg / m^3]
rho_Gl = 3860 ;
cp_Ti = 523;            %Specific heat capacity [J/(kg K)]
cp_Gl = 670;
k_Ti = 17 ;              %Heat conductivity [W/(m K)]
k_Gl = 0.8 ;
alpha_newt = 100;   %[W/ (m2 K)]    Newton convection

thickness = 1/100;      %Thickness (0.01 m)

Tc = 20;        %Temperature control rightmost boundary

TinfD = 40;     %Day time ambient air on leftmost boundary

TinfN = -96;    %Night time ambient air on leftmost boundary

coord = p'/100;         % coordinates of nodes [m]
enod=t(1:3,:)';         % nodes of elements
nelm=size(enod,1);      % number of elements
nen = 3;                % number of nodes per element
nnod=size(coord,1);     % number of nodes
dof=(1:nnod)';          % give each dof a number, dof for heat problem
ndof = length(dof);     % 1 degree of freedom per node
dof_S=[(1:nnod)',(nnod+1:2*nnod)']; %dof_S for mechanical problem
edof = zeros(nelm, 4);  %Preallocate memory
edof_S = zeros(nelm, 7);

for ie=1:nelm
    edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];     %For the mechanical problem
    edof(ie,:)=[ie,enod(ie,:)];                                                         %edof for heat equation
end

% Check which segments that should have convections
er = e([1 2 5],:);                  % Reduced e
conv_segments = [4 23 31 11];       % Choosen boundary segments
edges_conv = [];
for i = 1:size(er,2)
    if ismember(er(3,i),conv_segments)
        edges_conv = [edges_conv er(1:2,i)];
    end
end

[ex, ey] = coordxtr(edof, coord, dof, nen);

%% K, stiffness for heat problem and C for transient

D_Ti = k_Ti*eye(2);     %consitutive matrix for titanium
D_Gl = k_Gl*eye(2);     %consitutive matrix for glas
K = zeros(ndof);        %preallocate memory
C = zeros(ndof);        %preallocate memory
ep = thickness;

for elnr = 1:nelm
    if t(4, elnr) == 1      %Check if element is titanium
        Ke = flw2te(ex(elnr,:), ey(elnr,:), ep, D_Ti);  %Create element stiffness matrix
        Ce = plantml(ex(elnr,:), ey(elnr,:), cp_Ti * rho_Ti*ep);    %create elemenet matrix for C in transient
    else                    %Has to be glass else
        Ke = flw2te(ex(elnr,:), ey(elnr,:), ep, D_Gl);
        Ce = plantml(ex(elnr,:), ey(elnr,:), cp_Gl * rho_Gl*ep);
    end
    indx = edof(elnr,2:end);    
    K(indx,indx) = K(indx,indx) + Ke;   %Assemble K
    C(indx, indx) = C(indx, indx) + Ce; %assemble C
end 

% Daytime stationary, stemp2te adds convection
[K_D, a_D, f_D] = stemp2te(K, edges_conv, p, ep, alpha_newt, dof, edof, TinfD, Tc);
caxis([-60 35])
% Night time, add convection
[K_N, a_N, f_N] = stemp2te(K, edges_conv, p, ep, alpha_newt, dof, edof, TinfN, Tc);
caxis([-60 35])

%a_D has temperature distribution during day. a_N during night
maxDayTemp = max(a_D)
maxNightTemp = max(a_N)

%% Transient heat
%Transient heat, start with night time stationary distribution
%but apply daytime convection

%tstep(C, K_D, a_N, f_D, edof, ex, ey);

%Transient heat, start with daytime stationary distribution
%but apply night time convection
hold on
tstep(C, K_N, a_D, f_N, edof, ex, ey);


%% c) Geometry stationary and transient and mechanical

%plante.m   good functions to know
%plantf.m
%plants.m

%Ts = a_D;  %temperatures during day
Ts = a_N;  %temperatures during night
mag = 100; % Magnification (due to small deformations)
%mag = 20;  %mag during night
T0 = 20;

coord = p'/100;
enod=t(1:3,:)'; % nodes of elements
nelm=size(enod,1); % number of elements
nnod=size(coord,1); % number of nodes
dof=(1:nnod)'; % dof number is node number
dof_S=[(1:nnod)',(nnod+1:2*nnod)']; % give each dof a number
ndof = length(dof_S);   %Is actually half of number of degrees of freedom
edof = zeros(nelm, 4);  %Preallocate memory
edof_S = zeros(nelm, 7);%Preallocate memory
ep = [2 thickness];     %t is thickness and 2 is specification that it is a strain problem

for ie=1:nelm   %Create edof matrices, edof_S is for mechanical problem
    edof_S(ie,:)=[ie dof_S(enod(ie,1),:), dof_S(enod(ie,2),:),dof_S(enod(ie,3),:)];
    edof(ie,:)=[ie,enod(ie,:)];
end

[ex, ey] = coordxtr(edof, coord, dof, nen);
D_Ti = isoDte(ny_Ti, E_Ti);     %Creates the D matrix for titanium/glas, 3x3 matrix
D_Gl = isoDte(ny_Gl, E_Gl);
eT = extract(edof, Ts);
K = zeros(2*ndof);
f = zeros(2*ndof, 1);

%Create stiffness matrix K and force vector f
for elnr = 1:nelm
    dT = mean(eT(elnr,:)) - T0;     %Mean of temperature of nodes in triangle element, stress free at 20
    if t(4, elnr) == 1  %Titanium part
        epsilon0 = (1+ny_Ti) * alpha_Ti * dT * [1 1 0]';
        Ke = plante(ex(elnr,:), ey(elnr,:), ep, D_Ti);                          %Stress in z
        thermal = plantf(ex(elnr,:), ey(elnr,:), ep, (D_Ti*epsilon0)');
    else    %Glass otherwise
        epsilon0 =  (1+ny_Ti) * alpha_Gl * dT * [1 1 0]';
        Ke = plante(ex(elnr,:), ey(elnr,:), ep, D_Gl);              
        thermal = plantf(ex(elnr,:), ey(elnr,:), ep, (D_Gl*epsilon0)' );
    end

    indx = edof_S(elnr,2:end);
    K(indx,indx) = K(indx,indx) + Ke;
    f(indx) = f(indx) + thermal;
end

% Check which segments that should have fixed displacement in x or y
er = e([1 2 5],:); % Reduced e
fixed_segments_y = [20 21 3 15 16 17 18 19]; % Choosen boundary segments
fixed_segments_x = [11]; 
bc = [];            %Medvetna om att denna lösningen generar dubbletter
for i = 1:size(er,2)
    if ismember(er(3,i),fixed_segments_y)   %Checks if segment is part of boundary segments
        bc = [bc [er(1:2,i)'+ndof ; 0 0]];  %
    elseif ismember(er(3,i),fixed_segments_x)
        bc = [bc [er(1:2,i)' ; 0 0]];    
    end
end

a = solveq(K, f, bc');       %Nodal displacement

eD = extract(edof_S, a);    %pairs nodes and displacement with elements

[ex, ey] = coordxtr(edof_S, coord, dof_S, nen);
% Calculate displaced coordinates
exd = ex + mag*eD(:,1:2:end);
eyd = ey + mag*eD(:,2:2:end);
figure()
hold on
patch(ex',ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
patch(ex',-ey',[0 0 0],'EdgeColor','none','FaceAlpha',0.3)
patch(exd',eyd',[0 0 0],'FaceAlpha',0.3)
patch(exd',-eyd',[0 0 0],'FaceAlpha',0.3)
axis equal
xlabel('x-position [m]')
ylabel('y-position [m]')
%title("Displacement field during day conditions[Magnitude enhancement " + mag +  "]")
title("Displacement field during night conditions[Magnitude enhancement " + mag +  "]")

%preallocate memoey for von Mises stress in element
Seff_el = zeros(nelm, 1);

for elnr = 1:nelm
    dT = mean(eT(elnr,:)) - T0;
    if t(4, elnr) == 1
        epsilon0 = alpha_Ti * dT * [1 1 0]';
        [es, et] = plants(ex(elnr,:), ey(elnr,:), ep, D_Ti, eD(elnr, :));
        es = es - (D_Ti*epsilon0)';
        sigz = ny_Ti*(es(1)+es(2)) - alpha_Ti*E_Ti*dT;
    else 
        epsilon0 = alpha_Gl * dT * [1 1 0]';
        [es, et] = plants(ex(elnr,:), ey(elnr,:), ep, D_Gl, eD(elnr, :));
        es = es - (D_Gl*epsilon0)';
        sigz = ny_Gl*(es(1)+es(2)) - alpha_Gl*E_Gl*dT;
    end
    %indx = edof_S(elnr,2:end);
    sigx = es(1);
    sigy = es(2);
    gamxy = es(3);
    Seff_el(elnr) = Seff_el(elnr) + sqrt( sigx^2 + sigy^2 + sigz^2 - sigx*sigy - sigx*sigz - sigy*sigz + 3*gamxy^2);
end

Seff_nod = zeros(ndof);
for i = 1:size(coord,1)
     [c0,c1] = find(edof(:,2:4)==i);
     Seff_nod(i,1) = sum(Seff_el(c0))/size(c0,1);    %Seff_nod = von Mises effective stress at nodes, Seff_el = --- in elements
end

[ex, ey] = coordxtr(edof, coord, dof, nen);
eT = extract(edof, Seff_nod);
figure()
hold on
patch(ex',ey',eT', 'EdgeColor','none')
patch(ex',-ey',eT', 'EdgeColor','none')
xlabel('x-position [m]')
ylabel('y-position [m]')
title('von Mises stress in camera lens during day conditions [Pa]');
%title('von Mises stress in camera lens during night conditions [Pa]');
colormap(jet);
colorbar;

%% d)
%a is displacement vector

%calculate T acording to \int N^TN dA and find magnitude of 
% \int u^T*u dA which is square sum of displacement
T = zeros(2*ndof);
for elnr = 1:nelm
    if t(4, elnr) == 1      %Check if element is titanium
        Te = plantml(ex(elnr,:), ey(elnr,:), thickness);
    else                    %Has to be glass else
        Te = plantml(ex(elnr,:), ey(elnr,:), thickness);
    end
    indx = edof(elnr,2:end);    
    T(indx, indx) = T(indx, indx) + Te;
end 

%dispMagD = a'*T*a;
dispMagN = a'*T*a;

%dispMagD
dispMagN


