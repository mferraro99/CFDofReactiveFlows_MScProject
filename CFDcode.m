% Project of Computational Fluid Dynamics of Reactive Flows, ay. 2021/2022
% Authors: 
% - Elisa Giovanna Faggioli, id. 10628955, 
%   elisagiovanna.faggioli@mail.polimi.it
% - Marcello Ferraro, id. 10602732, 
%   marcello.ferraro@mail.polimi.it
%
% TITLE: STUDY OF THE FLUID DYNAMICS AND REACTION IN A TUBULAR REACTOR

clc;
clear variables;
close all;

%% PLEASE, SELECT A CASE: 
% CASE 1: uflow top and bottom = 1 m/s (only species A) 
%         uflow central = 0.5 m/s (only species D)
% CASE 2: uflow top and bottom = 0.5 m/s (only species D) 
%         uflow central = 1 m/s (only species A)
%//////////////////////////////////////////////////////////////////////////
                               CASE = 1;
%//////////////////////////////////////////////////////////////////////////
if CASE == 1
    uflow_sidepipes = 1;     % Inlet flow velocity u, side pipes [m/s]
    uflow_centralpipe = 0.5; % Inlet flow velocity u, central pipe [m/s]
    xA_sidepipes = 1;        % Molar fraction (A) side pipes
    xA_centralpipe = 0;      % Molar fraction (A) central pipe
    xD_sidepipes = 0;        % Molar fraction (D) side pipes
    xD_centralpipe = 1;      % Molar fraction (D) central pipe
    beta = 1.7;              % SOR coefficient
    max_error = 1e-6;        % Error for convergence
end
if CASE == 2
    uflow_sidepipes = 0.5;   % Inlet flow velocity u, side pipes [m/s]
    uflow_centralpipe = 1;   % Inlet flow velocity u, central pipe [m/s]
    xA_sidepipes = 0;        % Molar fraction (A) side pipes
    xA_centralpipe = 1;      % Molar fraction (A) central pipe
    xD_sidepipes = 1;        % Molar fraction (D) side pipes
    xD_centralpipe = 0;      % Molar fraction (D) central pipe
    beta = 1.9;              % SOR coefficient
    max_error = 1e-6;        % Error for convergence
end
if CASE ~=1 && CASE ~=2
    error('ERROR: "CASE" IS NOT EQUAL TO 1 OR 2, CHOOSE A VALID SCENARIO');
end

%% Data
%---------------------------Computational data-----------------------------
Lx = 1;           % Length of the top pipe [m]
Ly = 0.1;         % Diameter of the pipe [m]
nx = 600;         % Number of physical points in x direction 
ny = 60;          % Number of physical points in y direction
hx = Lx/nx;       % Step in x direction [m]
hy = Ly/ny;       % Step in y direction [m]
tau = 4;          % Simulation time [s]
max_iter = 10000; % Max number of iterations (Poisson solution)

%-----------------Fluidodynamic data (for all species)---------------------
MW = 46;              % Molecular weigth [kg/kmol]
T = 70+273.15;        % Temperature of the system [K]
P = 101325;           % Pressure of the system [Pa]
R = 8.314;            % Constant of ideal gases [J/mol/K]
rho = P*MW/R/T*10^-3; % Density of the system [kg/m3]
mu = 1e-4;            % Dynamic viscosity [Pa s]
nu = mu/rho;          % Kinematic viscosity [m2/s]
Diff = 4e-5;          % Diffusivity [m2/s]

%-----------------------------Kinetic data---------------------------------
k1 = 3.5;        % Kinetic constant of reaction 1 [1/s]
k2 = 1.2;        % Kinetic constant of reaction 2 [1/s]
k3 = 95e-6;      % Kinetic constant of reaction 3 [m3/mol/s]
Ctot0 = P/R/T;   % Total inlet concentration in every pipe [mol/m3]
nuA = [-1 0 0];  % Stoichiometric matrix (A)
nuB = [1 -1 -1]; % Stoichiometric matrix (B)
nuC = [0 1 0];   % Stoichiometric matrix (C)
nuD = [0 0 -1];  % Stoichiometric matrix (D)
nuE = [0 0 2];   % Stoichiometric matrix (E)

%% Preliminary operations

%-------Wall velocities (for reflection method boundary conditions)--------
un = 0; % North wall velocity [m/s]
us = 0; % South wall velocity [m/s]
ve = 0; % East wall velocity [m/s]
vw = 0; % West wall velocity [m/s]

%------------------------Time step calculation-----------------------------
sigma = 0.4; % Safety factor (to ensure stability)
% Since we are calculating the time increment at inlet velocity, we must
% add a safety factor to take into account the increment in velocity inside
% the reactor (which is unknown a priori), unless we use an adaptive time
% step increment

%--------------Time step rquired by Navier-Stokes stability----------------  
dt_conv_ns = min(4*nu/uflow_sidepipes^2,4*nu/uflow_centralpipe^2);                  
dt_diff_ns = min(hx^2/4/nu,hy^2/4/nu);   
dt_ns = min(dt_diff_ns,dt_conv_ns);

%-----------Time step required by transport of species stability-----------
dt_diff_sp = min(hx^2/4/Diff,hy^2/4/Diff);
dt_conv_sp = min(4*Diff/uflow_sidepipes^2,4*Diff/uflow_centralpipe^2);
dt_sp = min(dt_diff_sp,dt_conv_sp);

%-------------------------Time step for transport--------------------------
dt = sigma*min(dt_ns,dt_sp);
nsteps = tau/dt;

%---Time step required by reaction terms (because of operator splitting)---
dt_react = min(1/k1 , min(1/k2,1/k3/2/Ctot0) );
nsteps_react = ceil(dt/dt_react);
dt_react = dt/(nsteps_react);

%----------Check for grid consistency with the reactor geometry------------
if (mod(nx,10)~=0 || nx<400) %In order to have grid insensitivity
    error('Requirements for nx : minimun 400, multiple of 10');
end

if (mod(ny,4)~=0 || mod(ny,10)~=0 || ny<40) %In order to have grid insensitivity
    error('Requirements for ny : minimum 40, multiple of 4 and 10');
end
%//////////////////////////////////////////////////////////////////////////
%                           SOLUTION ALGORITHM                   
%//////////////////////////////////////////////////////////////////////////

%% MEMORY ALLOCATION: MAIN FIELD
u = zeros(nx+1,ny+2);
v = zeros(nx+2,ny+1);
p = zeros(nx+2,ny+2);
CA = zeros(nx+2,ny+2);
CB = zeros(nx+2,ny+2);
CC = zeros(nx+2,ny+2);
CD = zeros(nx+2,ny+2);
CE = zeros(nx+2,ny+2);

%% MEMORY ALLOCATION: TOP AND BOTTOM PIPE (COMPUTING ONLY ONE TIME)
% Subscript (s): side
%----------------------------Grid definition-------------------------------
nx_s = 2/10*nx; % Number of x points top and bottom pipe
ny_s = 1/10*ny; % Number of y points top and bottom pipe

%----------------------------Main field side-------------------------------
u_side = zeros(nx_s+1,ny_s+2);
v_side = zeros(nx_s+2,ny_s+1);
p_side = zeros(nx_s+2,ny_s+2);
CA_side = zeros(nx_s+2,ny_s+2);
CB_side = zeros(nx_s+2,ny_s+2);
CC_side = zeros(nx_s+2,ny_s+2);
CD_side = zeros(nx_s+2,ny_s+2);
CE_side = zeros(nx_s+2,ny_s+2);
CAstar_side =zeros(nx_s+2,ny_s+2);
CBstar_side = zeros(nx_s+2,ny_s+2);
CCstar_side = zeros(nx_s+2,ny_s+2);
CDstar_side = zeros(nx_s+2,ny_s+2);
CEstar_side = zeros(nx_s+2,ny_s+2);

%-----------------------Temporary velocity field side----------------------
ut_side=zeros(nx_s+1,ny_s+2);
vt_side=zeros(nx_s+2,ny_s+1);

%-------------------Gamma coefficients for Poisson equation----------------
gammax_side = zeros(nx_s+2,ny_s+2) + 2;
gammax_side(2:nx_s+1,2) = 1;            % South wall
gammax_side(2:nx_s+1,ny_s+1) = 1;       % North wall 

gammay_side = zeros(nx_s+2,ny_s+2)+2; 
gammay_side(2,2:ny_s+1) = 1;            % West boundary

%% MEMORY ALLOCATION: CENTRAL PIPE
% Subscript (c): central
%----------------------------Grid definition-------------------------------
nx_c = 2/10*nx; % Number of x points central pipe
ny_c = 1/2*ny;  % Number of y points central pipe

%---------------------------Main field central-----------------------------
u_central = zeros(nx_c+1,ny_c+2);
v_central = zeros(nx_c+2,ny_c+1);
p_central = zeros(nx_c+2,ny_c+2);
CA_central = zeros(nx_c+2,ny_c+2);
CB_central = zeros(nx_c+2,ny_c+2);
CC_central = zeros(nx_c+2,ny_c+2);
CD_central = zeros(nx_c+2,ny_c+2);
CE_central = zeros(nx_c+2,ny_c+2);
CAstar_central = zeros(nx_c+2,ny_c+2);
CBstar_central = zeros(nx_c+2,ny_c+2);
CCstar_central = zeros(nx_c+2,ny_c+2);
CDstar_central = zeros(nx_c+2,ny_c+2);
CEstar_central = zeros(nx_c+2,ny_c+2);

%---------------------Temporary velocity field central---------------------
ut_central=zeros(nx_c+1,ny_c+2);
vt_central=zeros(nx_c+2,ny_c+1);

%-------------------Gamma coefficients for Poisson equation----------------
gammax_central = zeros(nx_c+2,ny_c+2) + 2; 
gammax_central(2:nx_c+1,2) = 1;            % South wall
gammax_central(2:nx_c+1,ny_c+1) = 1;       % North wall 

gammay_central = zeros(nx_c+2,ny_c+2)+2; 
gammay_central(2,2:ny_c+1) = 1;            % West boundary

%% MEMORY ALLOCATION: MIXING CHAMBER
% Subscript (b): box
%----------------------------Grid definition-------------------------------
nx_b = 2/10*nx; % Number of x points mixing chamber
ny_b = ny;      % Number of y points mixing chamber

%------------------------Main field mixing chamber-------------------------
u_box = zeros(nx_b+1,ny_b+2);
v_box = zeros(nx_b+2,ny_b+1);
p_box = zeros(nx_b+2,ny_b+2);
CA_box = zeros(nx_b+2,ny_b+2);
CB_box = zeros(nx_b+2,ny_b+2);
CC_box = zeros(nx_b+2,ny_b+2);
CD_box = zeros(nx_b+2,ny_b+2);
CE_box = zeros(nx_b+2,ny_b+2);
CAstar_box = zeros(nx_b+2,ny_b+2);
CBstar_box = zeros(nx_b+2,ny_b+2);
CCstar_box = zeros(nx_b+2,ny_b+2);
CDstar_box = zeros(nx_b+2,ny_b+2);
CEstar_box = zeros(nx_b+2,ny_b+2);

%------------------Temporary velocity field mixing chamber-----------------
ut_box = zeros(nx_b+1,ny_b+2);
vt_box = zeros(nx_b+2,ny_b+1);

%-------------------Gamma coefficients for Poisson equation----------------
gammax_box = zeros(nx_b+2,ny_b+2) + 2; 
gammax_box(2:nx_b+1,2) = 1;               % South wall 
gammax_box(2:nx_b+1,ny_b+1) = 1;          % North wall 

gammay_box = zeros(nx_b+2,ny_b+2)+2; 
gammay_box(2,2:ny_b+1) = 1;               % West wall/boundary
gammay_box(nx_b+1,3/4*ny_b+2:ny_b+1) = 1; % East wall, over outlet pipe 
gammay_box(nx_b+1,2:ny_b/4+1) = 1;        % East wall, under outlet pipe

%% MEMORY ALLOCATION: OUTLET PIPE
%----------------------------Grid definition-------------------------------
nx_o = 6/10*nx; % Number of x points outlet pipe
ny_o = 1/2*ny;  % Number of y points outlet pipe

%-------------------------Main field outlet pipe---------------------------
u_out = zeros(nx_o+1,ny_o+2);
v_out = zeros(nx_o+2,ny_o+1);
p_out = zeros(nx_o+2,ny_o+2);
CA_out = zeros(nx_o+2,ny_o+2);
CB_out = zeros(nx_o+2,ny_o+2);
CC_out = zeros(nx_o+2,ny_o+2);
CD_out = zeros(nx_o+2,ny_o+2);
CE_out = zeros(nx_o+2,ny_o+2);
CAstar_out = zeros(nx_o+2,ny_o+2);
CBstar_out = zeros(nx_o+2,ny_o+2);
CCstar_out = zeros(nx_o+2,ny_o+2);
CDstar_out = zeros(nx_o+2,ny_o+2);
CEstar_out = zeros(nx_o+2,ny_o+2);

%-------------------Temporary velocity field outlet pipe-------------------
ut_out = zeros(nx_o+1,ny_o+2);
vt_out = zeros(nx_o+2,ny_o+1);

%-------------------Gamma coefficients for Poisson equation----------------
gammax_out = zeros(nx_o+2,ny_o+2)+2;
gammax_out(2:nx_o+1,2) = 1;          % South wall 
gammax_out(2:nx_o+1,ny_o+1) = 1;     % North wall

gammay_out = zeros(nx_o+2,ny_o+2)+2;
gammay_out(2,2:ny_o+1) = 1;          % West boundary


%% SOLUTION TIME LOOP
t = 0.0; %Set time equal to 0
for i = 1:nsteps

    %-------------------TOP AND BOTTOM PIPE SOLUTION-----------------------
    % 1) Projection algorithm

    % Boundary conditions
    u_side(1:nx_s+1,1) = 2*us - u_side(1:nx_s+1,2); % Reflection South 
    u_side(1:nx_s+1,ny_s+2)=2*un - u_side(1:nx_s+1,ny_s+1); % Reflection North
    v_side(1,1:ny_s+1) = 2*vw - v_side(2,1:ny_s+1); % Reflection West
    v_side(nx_s+2,1:ny_s+1) = 2*ve - v(nx_s+1,1:ny_s+1); % Reflection East
    u_side(1,2:ny_s+1) = uflow_sidepipes; % Dirichlet inlet flow on u
    u_side(nx_s+1,2:ny_s+1) = u_side(nx_s,2:ny_s+1); % Neumann outlet flow on u
    v_side(nx_s+2,2:ny_s+1) = v_side(nx_s+1,2:ny_s+1); % Neumann outlet flow on v

    % Advection-Diffusion equation
    [ut_side , vt_side] = AdvectionDiffusion2D(ut_side, vt_side, u_side,...
                                 v_side, nx_s, ny_s, hx, hy, dt, nu, Diff);

    % Assign inlet/outlet velocities at the temporary fields
    ut_side(1,2:ny_s+1) = uflow_sidepipes;
    ut_side(nx_s+1,2:ny_s+1) = u_side(nx_s+1,2:ny_s+1);
    vt_side(nx_s+2,2:ny_s+1) = v_side(nx_s+2,2:ny_s+1);

    % Pressure equation (Poisson)
    [p_side , iter_side] = Poisson2D(p_side, ut_side, vt_side, nx_s,...
      ny_s, hx, hy, dt, beta, max_iter, max_error,gammax_side,gammay_side);
                           
    % Velocity correction
    u_side(2:nx_s,2:ny_s+1)=ut_side(2:nx_s,2:ny_s+1)-...
               (dt/hx)*(p_side(3:nx_s+1,2:ny_s+1)-p_side(2:nx_s,2:ny_s+1));
    v_side(2:nx_s+1,2:ny_s)=vt_side(2:nx_s+1,2:ny_s)-...
               (dt/hy)*(p_side(2:nx_s+1,3:ny_s+1)-p_side(2:nx_s+1,2:ny_s));

    % 2) Transport of species

    % Boundary conditions
    CA_side = WallsBoundaryConditions(CA_side,nx_s,ny_s);
    CB_side = WallsBoundaryConditions(CB_side,nx_s,ny_s);
    CC_side = WallsBoundaryConditions(CC_side,nx_s,ny_s);
    CD_side = WallsBoundaryConditions(CD_side,nx_s,ny_s);
    CE_side = WallsBoundaryConditions(CE_side,nx_s,ny_s);
    CA_side(1,2:ny_s+1) = Ctot0*xA_sidepipes;
    CD_side(1,2:ny_s+1) = Ctot0*xD_sidepipes;
    CB_side(1,2:ny_s+1) = 0;
    CC_side(1,2:ny_s+1) = 0;
    CE_side(1,2:ny_s+1) = 0;

    % Advection-Diffusion for species
    CAstar_side=TemporaryConcentrations(CA_side,u_side,v_side,dt,hx,hy,Diff,nx_s,ny_s);
    CBstar_side=TemporaryConcentrations(CB_side,u_side,v_side,dt,hx,hy,Diff,nx_s,ny_s);
    CCstar_side=TemporaryConcentrations(CC_side,u_side,v_side,dt,hx,hy,Diff,nx_s,ny_s);
    CDstar_side=TemporaryConcentrations(CD_side,u_side,v_side,dt,hx,hy,Diff,nx_s,ny_s);
    CEstar_side=TemporaryConcentrations(CE_side,u_side,v_side,dt,hx,hy,Diff,nx_s,ny_s);

    % 3) Reaction of species

    % Calculate the rate of reaction (local) in every grid point
    r1 = CAstar_side*k1;
    r2 = CBstar_side*k2;
    r3 = CBstar_side.*CDstar_side*k3;

    % Reaction step
    CA_side = BatchSolution(CAstar_side,nx_s,ny_s,dt_react,nuA,r1,r2,r3);
    CB_side = BatchSolution(CBstar_side,nx_s,ny_s,dt_react,nuB,r1,r2,r3);
    CC_side = BatchSolution(CCstar_side,nx_s,ny_s,dt_react,nuC,r1,r2,r3);
    CD_side = BatchSolution(CDstar_side,nx_s,ny_s,dt_react,nuD,r1,r2,r3);
    CE_side = BatchSolution(CEstar_side,nx_s,ny_s,dt_react,nuE,r1,r2,r3);

    %------------------------CENTRAL PIPE SOLUTION-------------------------
    % 1) Projection algorithm

    % Boundary conditions
    u_central(2:nx_c+1,1) = 2*us - u_central(2:nx_c+1,2); % Reflection South
    u_central(2:nx_c+1,ny_c+2)=2*un - u_central(2:nx_c+1,ny_c+1); % Reflection North
    v_central(1:1:ny_c+1) = 2*vw - v_central(2,1:ny_c+1); % Reflection West
    v_central(nx_c+2,1:ny_c+1) = 2*ve - v_central(nx_c+1,1:ny_c+1); % Reflection East
    u_central(1,2:ny_c+1) = uflow_centralpipe; % Dirichlet inlet flow
    u_central(nx_c+1,2:ny_c+1) = u_central(nx_c,2:ny_c+1); % Neumann outlet flow on u
    v_central(nx_c+2,2:ny_c+1) = v_central(nx_c+1,2:ny_c+1); % Neumann outlet flow on v

    % Advection-Diffusion equation
    [ut_central , vt_central] = AdvectionDiffusion2D(ut_central, vt_central,...
                   u_central, v_central, nx_c, ny_c, hx, hy, dt, nu, Diff);

    % Assign inlet/outlet velocities at the temporary fields
    ut_central(1,2:ny_c+1) = uflow_centralpipe;
    ut_central(nx_c+1,2:ny_c+1) = u_central(nx_c+1,2:ny_c+1);
    vt_central(nx_c+2,2:ny_c+1) = v_central(nx_c+2,2:ny_c+1);

    % Pressure equation (Poisson)
    [p_central , iter_central] = Poisson2D(p_central, ut_central, vt_central,...
    nx_c, ny_c, hx, hy, dt, beta, max_iter, max_error,gammax_central,gammay_central);

    % Velocity correction
    u_central(2:nx_c,2:ny_c+1)=ut_central(2:nx_c,2:ny_c+1)-...
        (dt/hx)*(p_central(3:nx_c+1,2:ny_c+1)-p_central(2:nx_c,2:ny_c+1));
    v_central(2:nx_c+1,2:ny_c)=vt_central(2:nx_c+1,2:ny_c)-...
        (dt/hy)*(p_central(2:nx_c+1,3:ny_c+1)-p_central(2:nx_c+1,2:ny_c));
                     
    % 2) Transport of species

    % Boundary conditions
    CA_central = WallsBoundaryConditions(CA_central,nx_c,ny_c);
    CB_central = WallsBoundaryConditions(CB_central,nx_c,ny_c);
    CC_central = WallsBoundaryConditions(CC_central,nx_c,ny_c);
    CD_central = WallsBoundaryConditions(CD_central,nx_c,ny_c);
    CE_central = WallsBoundaryConditions(CE_central,nx_c,ny_c);
    CA_central(1,1:ny_c+1) = Ctot0*xA_centralpipe;
    CD_central(1,1:ny_c+1) = Ctot0*xD_centralpipe;
    CB_central(1,1:ny_c+1) = 0;
    CC_central(1,1:ny_c+1) = 0;
    CE_central(1,1:ny_c+1) = 0;

    % Advection-Diffusion for species
    CAstar_central=TemporaryConcentrations(CA_central,u_central,v_central,dt,hx,hy,Diff,nx_c,ny_c);
    CBstar_central=TemporaryConcentrations(CB_central,u_central,v_central,dt,hx,hy,Diff,nx_c,ny_c);
    CCstar_central=TemporaryConcentrations(CC_central,u_central,v_central,dt,hx,hy,Diff,nx_c,ny_c);
    CDstar_central=TemporaryConcentrations(CD_central,u_central,v_central,dt,hx,hy,Diff,nx_c,ny_c);
    CEstar_central=TemporaryConcentrations(CE_central,u_central,v_central,dt,hx,hy,Diff,nx_c,ny_c);

    % 3) Reaction of species

    % Calculate the rate of reaction (local) in every grid point
    r1 = CAstar_central*k1;
    r2 = CBstar_central*k2;
    r3 = CBstar_central.*CDstar_central*k3;

    % Reaction step
    CA_central = BatchSolution(CAstar_central,nx_c,ny_c,dt_react,nuA,r1,r2,r3);
    CB_central = BatchSolution(CBstar_central,nx_c,ny_c,dt_react,nuB,r1,r2,r3);
    CC_central = BatchSolution(CCstar_central,nx_c,ny_c,dt_react,nuC,r1,r2,r3);
    CD_central = BatchSolution(CDstar_central,nx_c,ny_c,dt_react,nuD,r1,r2,r3);
    CE_central = BatchSolution(CEstar_central,nx_c,ny_c,dt_react,nuE,r1,r2,r3);

    %------------------------MIXING CHAMBER SOLUTION-----------------------
    % 1) Projection algorithm

    % Boundary conditions
    u_box(1:nx_b+1,1) = 2*us - u_box(1:nx_b+1,2); % Reflection South
    u_box(1:nx_b+1,ny_b+2) = 2*un - u_box(1:nx_b+1,ny_b+1); % Reflection North
    v_box(1,1:ny_b+1) = 2*vw - v_box(2,1:ny_b+1); % Reflection West
    v_box(nx_b+2,1:ny_b+1) = 2*ve - v_box(nx_b+1,1:ny_b+1); % Reflection East

    u_box(1,2:ny_b/10+1) = u_side(nx_s+1,2:ny_s+1); % Dirichlet inlet u bottom pipe
    v_box(1,2:ny_b/10) = v_side(nx_s+1,2:ny_s); % Dirichlet inlet v bottom pipe

    u_box(1,ny_b/4+2:3/4*ny_b+1) = u_central(nx_c+1,2:ny_c+1); % Dirichlet inlet u central pipe
    v_box(1,ny_b/4+2:3/4*ny_b) = v_central(nx_c+1,2:ny_c); % Dirichlet inlet v central pipe

    u_box(1,9/10*ny_b+2:ny_b+1) = u_side(nx_s+1,2:ny_s+1); % Dirichlet inlet u top pipe
    v_box(1,9/10*ny_b+2:ny_b) = v_side(nx_s+1,2:ny_s); % Dirichlet inlet v top pipe

    u_box(nx_b+1,ny_b/4+2:3/4*ny_b+1) = u_box(nx_b,ny_b/4+2:3/4*ny_b+1); % Neumann outlet u
    v_box(nx_b+2,ny_b/4+2:3/4*ny_b) = v_box(nx_b+1,ny_b/4+2:3/4*ny_b); % Neumann outlet v

    % Advection-Diffusion equation
    [ut_box , vt_box] = AdvectionDiffusion2D(ut_box, vt_box, u_box, v_box,...
                                         nx_b, ny_b, hx, hy, dt, nu, Diff);

    % Assign inlet/outlet velocities at the temporary fields
    ut_box(1,2:ny_b/10+1) = u_side(nx_s+1,2:ny_s+1); 
    vt_box(1,2:ny_b/10) = v_side(nx_s+1,2:ny_s); 
    ut_box(1,ny_b/4+2:3/4*ny_b+1) = u_central(nx_c+1,2:ny_c+1); 
    vt_box(1,ny_b/4+2:3/4*ny_b) = v_central(nx_c+1,2:ny_c); 
    ut_box(1,9/10*ny_b+2:ny_b+1) = u_side(nx_s+1,2:ny_s+1); 
    vt_box(1,9/10*ny_b+2:ny_b) = v_side(nx_s+1,2:ny_s);    
    
    vt_box(nx_b+2,ny_b/4+2:3/4*ny_b) = v_box(nx_b+1,ny_b/4+2:3/4*ny_b); 
    ut_box(nx_b+1,ny_b/4+2:3/4*ny_b+1) = u_box(nx_b,ny_b/4+2:3/4*ny_b+1); 
    
    % Pressure equation (Poisson)
    [p_box , iter_box] = Poisson2D(p_box, ut_box, vt_box, nx_b, ny_b, hx, hy, dt, ...
                           beta, max_iter, max_error,gammax_box,gammay_box);

    % Velocity correction
    u_box(2:nx_b,2:ny_b+1)=ut_box(2:nx_b,2:ny_b+1)-...
                 (dt/hx)*(p_box(3:nx_b+1,2:ny_b+1)-p_box(2:nx_b,2:ny_b+1));
    v_box(2:nx_b+1,2:ny_b)=vt_box(2:nx_b+1,2:ny_b)-...
                 (dt/hy)*(p_box(2:nx_b+1,3:ny_b+1)-p_box(2:nx_b+1,2:ny_b));

    % 2) Transport of species

    % Boundary conditions
    CA_box = WallsBoundaryConditions(CA_box,nx_b,ny_b);
    CB_box = WallsBoundaryConditions(CB_box,nx_b,ny_b);
    CC_box = WallsBoundaryConditions(CC_box,nx_b,ny_b);
    CD_box = WallsBoundaryConditions(CD_box,nx_b,ny_b);
    CE_box = WallsBoundaryConditions(CE_box,nx_b,ny_b);
    
    CA_box(1,2:ny_b/10+1) = CA_side(nx_s+1,2:ny_s+1); %BOTTOM INLET
    CB_box(1,2:ny_b/10+1) = CB_side(nx_s+1,2:ny_s+1); 
    CC_box(1,2:ny_b/10+1) = CC_side(nx_s+1,2:ny_s+1); 
    CD_box(1,2:ny_b/10+1) = CD_side(nx_s+1,2:ny_s+1); 
    CE_box(1,2:ny_b/10+1) = CE_side(nx_s+1,2:ny_s+1); 

    CA_box(1,ny_b/4+2:3/4*ny_b+1) = CA_central(nx_c+1,2:ny_c+1); %CENTRAL INLET
    CB_box(1,ny_b/4+2:3/4*ny_b+1) = CB_central(nx_c+1,2:ny_c+1); 
    CC_box(1,ny_b/4+2:3/4*ny_b+1) = CC_central(nx_c+1,2:ny_c+1); 
    CD_box(1,ny_b/4+2:3/4*ny_b+1) = CD_central(nx_c+1,2:ny_c+1); 
    CE_box(1,ny_b/4+2:3/4*ny_b+1) = CE_central(nx_c+1,2:ny_c+1); 

    CA_box(1,9/10*ny_b+2:ny_b+1) = CA_side(nx_s+1,2:ny_s+1); %TOP INLET
    CB_box(1,9/10*ny_b+2:ny_b+1) = CB_side(nx_s+1,2:ny_s+1); 
    CC_box(1,9/10*ny_b+2:ny_b+1) = CC_side(nx_s+1,2:ny_s+1); 
    CD_box(1,9/10*ny_b+2:ny_b+1) = CD_side(nx_s+1,2:ny_s+1); 
    CE_box(1,9/10*ny_b+2:ny_b+1) = CE_side(nx_s+1,2:ny_s+1); 

    % Advection-Diffusion for species
    CAstar_box=TemporaryConcentrations(CA_box,u_box,v_box,dt,hx,hy,Diff,nx_b,ny_b);
    CBstar_box=TemporaryConcentrations(CB_box,u_box,v_box,dt,hx,hy,Diff,nx_b,ny_b);
    CCstar_box=TemporaryConcentrations(CC_box,u_box,v_box,dt,hx,hy,Diff,nx_b,ny_b);
    CDstar_box=TemporaryConcentrations(CD_box,u_box,v_box,dt,hx,hy,Diff,nx_b,ny_b);
    CEstar_box=TemporaryConcentrations(CE_box,u_box,v_box,dt,hx,hy,Diff,nx_b,ny_b);

    % 3) Reaction of species

    % Calculate the rate of reaction (local) in every grid point
    r1 = CAstar_box*k1;
    r2 = CBstar_box*k2;
    r3 = CBstar_box.*CDstar_box*k3;

    % Reaction step
    CA_box = BatchSolution(CAstar_box,nx_b,ny_b,dt_react,nuA,r1,r2,r3);
    CB_box = BatchSolution(CBstar_box,nx_b,ny_b,dt_react,nuB,r1,r2,r3);
    CC_box = BatchSolution(CCstar_box,nx_b,ny_b,dt_react,nuC,r1,r2,r3);
    CD_box = BatchSolution(CDstar_box,nx_b,ny_b,dt_react,nuD,r1,r2,r3);
    CE_box = BatchSolution(CEstar_box,nx_b,ny_b,dt_react,nuE,r1,r2,r3);

    %-------------------------OUTLET PIPE SOLUTION-------------------------
    % 1) Projection algorithm

    % Boundary conditions
    u_out(1:nx_o+1,1) = 2*us - u_out(1:nx_o+1,2); % Reflection South
    u_out(1:nx_o+1,ny_o+2) = 2*un - u_out(1:nx_o+1,ny_o+1); % Reflection North
    v_out(1,1:ny_o+1) = 2*vw - v_out(2,1:ny_o+1); % Reflection West
    v_out(nx_o+2,1:ny_o+1) = 2*ve - v_out(nx_o+1,1:ny_o+1); % Reflection East
    u_out(1,2:ny_o+1) = u_box(nx_b+1,ny_b/4+2:3*ny_b/4+1); % DIrichlet inlet flow u
    v_out(1,2:ny_o) = v_box(nx_b+1,ny_b/4+2:3*ny_b/4); % DIrichlet inlet flow v
    u_out(nx_o+1,2:ny_o+1) = u_out(nx_o,2:ny_o+1); % Neumann outlet flow u
    v_out(nx_o+2,2:ny_o) = v_out(nx_o+1,2:ny_o); % Neumann outlet flow v

    % Advection-Diffusion equation
    [ut_out , vt_out] = AdvectionDiffusion2D(ut_out, vt_out, u_out, v_out,...
                                         nx_o, ny_o, hx, hy, dt, nu, Diff);    

    % Assign inlet/outlet velocities at the temporary fields
    ut_out(1,2:ny_o+1) = u_box(nx_b+1,ny_b/4+2:3*ny_b/4+1); 
    vt_out(1,2:ny_o) = v_box(nx_b+1,ny_b/4+2:3*ny_b/4); 

    ut_out(nx_o+1,2:ny_o+1) = u_out(nx_o,2:ny_o+1); 
    vt_out(nx_o+2,2:ny_o) = v_out(nx_o+1,2:ny_o);
    
    % Pressure equation (Poisson)
    [p_out , iter_out] = Poisson2D(p_out, ut_out, vt_out, nx_o, ny_o, hx, hy, dt, ...
                           beta, max_iter, max_error,gammax_out,gammay_out);  

    % Velocity correction
    u_out(2:nx_o,2:ny_o+1)=ut_out(2:nx_o,2:ny_o+1)-...
                 (dt/hx)*(p_out(3:nx_o+1,2:ny_o+1)-p_out(2:nx_o,2:ny_o+1));
    v_out(2:nx_o+1,2:ny_o)=vt_out(2:nx_o+1,2:ny_o)-...
                 (dt/hy)*(p_out(2:nx_o+1,3:ny_o+1)-p_out(2:nx_o+1,2:ny_o));

    % 2) Transport of species

    % Boundary conditions
    CA_out = WallsBoundaryConditions(CA_out,nx_o,ny_o);
    CB_out = WallsBoundaryConditions(CB_out,nx_o,ny_o);
    CC_out = WallsBoundaryConditions(CC_out,nx_o,ny_o);
    CD_out = WallsBoundaryConditions(CD_out,nx_o,ny_o);
    CE_out = WallsBoundaryConditions(CE_out,nx_o,ny_o);

    CA_out(1,2:ny_o+1) = CA_box(nx_b+1,ny_b/4+2:3*ny_b/4+1); % INLET CONCENTRATIONS
    CB_out(1,2:ny_o+1) = CB_box(nx_b+1,ny_b/4+2:3*ny_b/4+1); 
    CC_out(1,2:ny_o+1) = CC_box(nx_b+1,ny_b/4+2:3*ny_b/4+1); 
    CD_out(1,2:ny_o+1) = CD_box(nx_b+1,ny_b/4+2:3*ny_b/4+1); 
    CE_out(1,2:ny_o+1) = CE_box(nx_b+1,ny_b/4+2:3*ny_b/4+1); 

    CAstar_out=TemporaryConcentrations(CA_out,u_out,v_out,dt,hx,hy,Diff,nx_o,ny_o);
    CBstar_out=TemporaryConcentrations(CB_out,u_out,v_out,dt,hx,hy,Diff,nx_o,ny_o);
    CCstar_out=TemporaryConcentrations(CC_out,u_out,v_out,dt,hx,hy,Diff,nx_o,ny_o);
    CDstar_out=TemporaryConcentrations(CD_out,u_out,v_out,dt,hx,hy,Diff,nx_o,ny_o);
    CEstar_out=TemporaryConcentrations(CE_out,u_out,v_out,dt,hx,hy,Diff,nx_o,ny_o);

    % 3) Reaction of species

    % Calculate the rate of reaction (local) in every grid point
    r1 = CAstar_out*k1;
    r2 = CBstar_out*k2;
    r3 = CBstar_out.*CDstar_out*k3;

    % Reaction step
    CA_out = BatchSolution(CAstar_out,nx_o,ny_o,dt_react,nuA,r1,r2,r3);
    CB_out = BatchSolution(CBstar_out,nx_o,ny_o,dt_react,nuB,r1,r2,r3);
    CC_out = BatchSolution(CCstar_out,nx_o,ny_o,dt_react,nuC,r1,r2,r3);
    CD_out = BatchSolution(CDstar_out,nx_o,ny_o,dt_react,nuD,r1,r2,r3);
    CE_out = BatchSolution(CEstar_out,nx_o,ny_o,dt_react,nuE,r1,r2,r3);

    %------------------------GRID RECONSTRUCTION---------------------------

    % Bottom pipe
    u(1:nx_s,1:ny/10+1) = u_side(1:nx_s,1:ny_s+1);
    v(1:nx_s+1,1:ny/10) = v_side(1:nx_s+1,1:ny_s);
    p(1:nx_s+1,1:ny/10+1) = p_side(1:nx_s+1,1:ny_s+1);

    CA(1:nx_s,1:ny/10+1) = CA_side(1:nx_s,1:ny_s+1);
    CB(1:nx_s,1:ny/10+1) = CB_side(1:nx_s,1:ny_s+1);
    CC(1:nx_s,1:ny/10+1) = CC_side(1:nx_s,1:ny_s+1);
    CD(1:nx_s,1:ny/10+1) = CD_side(1:nx_s,1:ny_s+1);
    CE(1:nx_s,1:ny/10+1) = CE_side(1:nx_s,1:ny_s+1);

    % Central pipe
    u(1:nx_s,ny/4+2:3/4*ny+1) = u_central(1:nx_s,2:ny_c+1);
    v(1:nx_s+1,ny/4+2:3/4*ny) = v_central(1:nx_s+1,2:ny_c);
    p(1:nx_s+1,ny/4+2:3/4*ny+1) = p_central(1:nx_s+1,2:ny_c+1);

    CA(1:nx_s,ny/4+2:3/4*ny+1) = CA_central(1:nx_s,2:ny_c+1);
    CB(1:nx_s,ny/4+2:3/4*ny+1) = CB_central(1:nx_s,2:ny_c+1);
    CC(1:nx_s,ny/4+2:3/4*ny+1) = CC_central(1:nx_s,2:ny_c+1);
    CD(1:nx_s,ny/4+2:3/4*ny+1) = CD_central(1:nx_s,2:ny_c+1);
    CE(1:nx_s,ny/4+2:3/4*ny+1) = CE_central(1:nx_s,2:ny_c+1);

    % Top pipe
    u(1:nx_s,9/10*ny+2:ny+2) = u_side(1:nx_s,2:ny_s+2);
    v(1:nx_s+1,9/10*ny+2:ny+1) = v_side(1:nx_s+1,2:ny_s+1);
    p(1:nx_s+1,9/10*ny+2:ny+2) = p_side(1:nx_s+1,2:ny_s+2);

    CA(1:nx_s,9/10*ny+2:ny+2) = CA_side(1:nx_s,2:ny_s+2);
    CB(1:nx_s,9/10*ny+2:ny+2) = CB_side(1:nx_s,2:ny_s+2);
    CC(1:nx_s,9/10*ny+2:ny+2) = CC_side(1:nx_s,2:ny_s+2);
    CD(1:nx_s,9/10*ny+2:ny+2) = CD_side(1:nx_s,2:ny_s+2);
    CE(1:nx_s,9/10*ny+2:ny+2) = CE_side(1:nx_s,2:ny_s+2);

    % Mixing chamber
    u(nx_s+1:nx_s+nx_b,:) = u_box(1:nx_b,:);
    v(nx_s+1:nx_s+nx_b+1,:) = v_box(1:nx_b+1,:);
    p(nx_s+2:nx_s+nx_b+1,:) = p_box(2:nx_b+1,:);

    CA(nx_s+1:nx_s+nx_b,:) = CA_box(1:nx_b,:);
    CB(nx_s+1:nx_s+nx_b,:) = CB_box(1:nx_b,:);
    CC(nx_s+1:nx_s+nx_b,:) = CC_box(1:nx_b,:);
    CD(nx_s+1:nx_s+nx_b,:) = CD_box(1:nx_b,:);
    CE(nx_s+1:nx_s+nx_b,:) = CE_box(1:nx_b,:);

    % Outlet pipe
    u(nx_s+nx_b+1:nx,ny/4+2:3/4*ny+1) = u_out(1:nx_o,2:ny_o+1);
    v(nx_s+nx_b+1:nx+1,ny/4+1:3/4*ny) = v_out(1:nx_o+1,2:ny_o+1);
    p(nx_s+nx_b+2:nx+1,ny/4+2:3/4*ny+1) = p_out(2:nx_o+1,2:ny_o+1);

    CA(nx_s+nx_b:nx+1,ny/4+2:3/4*ny+1) = CA_out(:,2:ny_o+1);
    CB(nx_s+nx_b:nx+1,ny/4+2:3/4*ny+1) = CB_out(:,2:ny_o+1);
    CC(nx_s+nx_b:nx+1,ny/4+2:3/4*ny+1) = CC_out(:,2:ny_o+1);
    CD(nx_s+nx_b:nx+1,ny/4+2:3/4*ny+1) = CD_out(:,2:ny_o+1);
    CE(nx_s+nx_b:nx+1,ny/4+2:3/4*ny+1) = CE_out(:,2:ny_o+1);

    %--------------------------TIME ADVANCING------------------------------
    
    if (mod(i,100)==1) %Print the number of iterations of the Poisson solution in every section
        fprintf('-----------------------------------------------------------------------------------------------------------');
        fprintf('\nStep: %d - Time: %f - Poisson iterations central - top/bottom - box - outlet: %d - %d - %d - %d\n', ...
                        i, t, iter_central, iter_side, iter_box, iter_out);
    end      
    t = t+dt; %Advance in time
end

%% POST PROCESSING
x=0:hx:Lx;
y=0:hy:Ly;
[X,Y]=meshgrid(x,y);        

% Field reconstruction
uu(1:nx+1,1:ny+1)=0.50*(u(1:nx+1,2:ny+2)+u(1:nx+1,1:ny+1));
vv(1:nx+1,1:ny+1)=0.50*(v(2:nx+2,1:ny+1)+v(1:nx+1,1:ny+1));
pp(1:nx+1,1:ny+1)=0.25*(p(1:nx+1,1:ny+1)+p(1:nx+1,2:ny+2)+...
                        p(2:nx+2,1:ny+1)+p(2:nx+2,2:ny+2));
ca(1:nx+1,1:ny+1)=0.25*(CA(1:nx+1,1:ny+1)+CA(1:nx+1,2:ny+2)+...
                        CA(2:nx+2,1:ny+1)+CA(2:nx+2,2:ny+2));
cb(1:nx+1,1:ny+1)=0.25*(CB(1:nx+1,1:ny+1)+CB(1:nx+1,2:ny+2)+...
                        CB(2:nx+2,1:ny+1)+CB(2:nx+2,2:ny+2));
cc(1:nx+1,1:ny+1)=0.25*(CC(1:nx+1,1:ny+1)+CC(1:nx+1,2:ny+2)+...
                        CC(2:nx+2,1:ny+1)+CC(2:nx+2,2:ny+2));
cd(1:nx+1,1:ny+1)=0.25*(CD(1:nx+1,1:ny+1)+CD(1:nx+1,2:ny+2)+...
                        CD(2:nx+2,1:ny+1)+CD(2:nx+2,2:ny+2));
ce(1:nx+1,1:ny+1)=0.25*(CE(1:nx+1,1:ny+1)+CE(1:nx+1,2:ny+2)+...
                        CE(2:nx+2,1:ny+1)+CE(2:nx+2,2:ny+2));

% Plot the maps for the concentrations
% Total concentration
subplot(231);
surface(X,Y,ca'+cb'+cc'+cd'+ce', 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
title('Total concentration [kmol/m3]','Interpreter','latex'); 
xlabel('x coordinate [m]','Interpreter','latex'); 
ylabel('y coordinate [m]','Interpreter','latex');
colorbar; shading interp;

% CA,CB,CC,CD,CE maps
subplot(232)
surface(X,Y,ca', 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
title('$C_{A}$ [kmol/m3]','Interpreter','Latex'); 
xlabel('x coordinate [m]','Interpreter','latex'); 
ylabel('y coordinate [m]','Interpreter','latex');
colorbar; shading interp;        

subplot(233)
surface(X,Y,cb', 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
title('$C_{B}$ [kmol/m3]','Interpreter','latex'); 
xlabel('x coordinate [m]','Interpreter','latex'); 
ylabel('y coordinate [m]','Interpreter','latex');
colorbar; shading interp;        

subplot(234)
surface(X,Y,cc', 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
title('$C_{C}$ [kmol/m3]','Interpreter','latex'); 
xlabel('x coordinate [m]','Interpreter','latex'); 
ylabel('y coordinate [m]','Interpreter','latex');
colorbar; shading interp;        

subplot(235)
surface(X,Y,cd', 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
title('$C_{D}$ [kmol/m3]','Interpreter','latex'); 
xlabel('x coordinate [m]','Interpreter','latex'); 
ylabel('y coordinate [m]','Interpreter','latex');
colorbar; shading interp;        

subplot(236)
surface(X,Y,ce', 'EdgeColor','none','LineStyle','none','FaceLighting','phong');
title('$C_{E}$ [kmol/m3]','Interpreter','latex'); 
xlabel('x coordinate [m]','Interpreter','latex'); 
ylabel('y coordinate [m]','Interpreter','latex');
colorbar; shading interp;

%% SPECIES AND PRODUCTION
maxCE=max(CE(end-1,:));
meanCE=mean(CE(end-1,:));

fprintf('\nThe maximum value for species E in case %d is %f\n',CASE,maxCE);
fprintf('The mean value for species E in case %d is %f\n',CASE,meanCE);

figure;
yy=linspace(0,0.1,length(CE(end,:)));
plot(CE(end-1,:),yy,'r','LineWidth',2)

% Determination of molar flux of E
% Flux in [mol/m^2/s]

Fe=(u(end-1,:)).*CE(end-1,:);
mean_FE=mean(Fe);

fprintf('The mean flux of E in case %d is %f\n',CASE,mean_FE)

prod=Fe*3600;
mean_P=mean(prod);

fprintf('The production of species E in case %d after 1h is %f\n',CASE,mean_P)
                            
%% FUNCTIONS
%//////////////////////////////////////////////////////////////////////////
%                               FUNCTIONS
%//////////////////////////////////////////////////////////////////////////

%------------------------POISSON EQUATION SOLVER---------------------------
function [p, iter] = Poisson2D( p, ut, vt, nx, ny, hx,hy, dt, ...
                            beta, max_iter, max_error,gammax,gammay)
    for iter=1:max_iter
        for i=2:nx+1
            for j=2:ny+1
                
                delta = (hy^2*( p(i+1,j)+p(i-1,j) )+ hx^2*( p(i,j+1)+p(i,j-1) ));
                S = -hx*hy^(2)/dt *( ut(i,j) - ut(i-1,j) ) - hx^(2)*hy/dt *( vt(i,j) - vt(i,j-1) );
                p(i,j) = beta/(gammay(i,j)*hy^2+gammax(i,j)*hx^2)*( delta+S )+(1-beta)*p(i,j);
                
            end
        end
    
        % Estimate the error
        epsilon=0.0; 
        for i=2:nx+1
            for j=2:ny+1
                delta = (hy^2*( p(i+1,j)+p(i-1,j) )+ hx^2*( p(i,j+1)+p(i,j-1) ));
                S = -hx*hy^(2)/dt *( ut(i,j) - ut(i-1,j) ) - hx^(2)*hy/dt *( vt(i,j) - vt(i,j-1) );             
                epsilon=epsilon+abs( p(i,j) - 1/(gammay(i,j)*hy^2+gammax(i,j)*hx^2)*( delta+S ) );
            end
        end
        epsilon = epsilon / (nx*ny);
    
        % Check the error
        if (epsilon <= max_error) % stop if converged
            break;
        end 
    end
end

%-----------------------ADVECTION-DIFFUSION EQUATION-----------------------
function [ut, vt] = AdvectionDiffusion2D( ut, vt, u, v, nx, ny, hx, hy, dt, nu ,Diff)

    % Temporary u-velocity
    for i=2:nx
        for j=2:ny+1 
            
            Pe_u = abs(u(i,j))*hx/Diff; %Check on u velocity for hybrid scheme
            if Pe_u>=2 %Upwind
                flux_e = 0.5*(u(i+1,j)   + u(i,j));
                flux_w = 0.5*(u(i,j)     + u(i-1,j));
                flux_n = 0.5*(v(i+1,j)   + v(i,j));
                flux_s = 0.5*(v(i+1,j-1) + v(i,j-1));

                if (flux_e > 0) 
                    ufe = u(i,j);   
                else 
                    ufe = u(i+1,j); 
                end
                if (flux_w > 0) 
                    ufw = u(i-1,j); 
                else 
                    ufw = u(i,j);   
                end
                if (flux_n > 0) 
                    ufn = u(i,j);   
                else 
                    ufn = u(i,j+1); 
                end
                if (flux_s > 0) 
                    ufs = u(i,j-1); 
                else 
                    ufs = u(i,j);   
                end

                ue = ufe^2;
                uw = ufw^2;
                un = flux_n*ufn;
                us = flux_s*ufs;
            else %Centered
                ue = 1/4 * ( u(i+1,j)+u(i,j) )^2;
                uw = 1/4 * ( u(i,j)+u(i-1,j) )^2;
                un = 1/4 * ( u(i,j+1)+u(i,j) )*( v(i+1,j)+v(i,j) );
                us = 1/4 * ( u(i,j)+u(i,j-1) )*( v(i+1,j-1)+v(i,j-1) );
            end
            
            A = (ue-uw)/hx +(un-us)/hy;
            D = nu/hx^2 * ( u(i+1,j)-2*u(i,j)+u(i-1,j) ) +...
                nu/hy^2 * ( u(i,j+1)-2*u(i,j)+u(i,j-1));
            
            ut(i,j)=u(i,j)+dt*(-A+D);   
        end
    end
    
    % Temporary v-velocity
    for i=2:nx+1
        for j=2:ny 
            Pe_v = abs(v(i,j))*hy/Diff; %Check on v velocity for hybrid scheme
            if Pe_v>=2 %Upwind
                flux_e = 0.5*(u(i,j)   + u(i,j+1));
                flux_w = 0.5*(u(i-1,j) + u(i-1,j+1));
                flux_n = 0.5*(v(i,j+1) + v(i,j));
                flux_s = 0.5*(v(i,j)   + v(i,j-1));

                if (flux_n > 0) 
                    vfn = v(i,j);   
                else 
                    vfn = v(i,j+1); 
                end
                if (flux_s > 0) 
                    vfs = v(i,j-1); 
                else 
                    vfs = v(i,j);   
                end
                if (flux_e > 0) 
                    vfe = v(i,j);   
                else 
                    vfe = v(i+1,j); 
                end
                if (flux_w > 0) 
                    vfw = v(i-1,j); 
                else 
                    vfw = v(i,j);   
                end

                vn = vfn^2;
                vs = vfs^2;
                ve = flux_e*vfe;
                vw = flux_w*vfw;
            else %Centered
                vn = 1/4 * ( v(i,j+1)+v(i,j) )^2;
                vs = 1/4 * ( v(i,j)+v(i,j-1) )^2;
                ve = 1/4 * ( u(i,j+1)+u(i,j) )*( v(i+1,j)+v(i,j) );
                vw = 1/4 * ( u(i-1,j+1)+u(i-1,j) )*( v(i,j)+v(i-1,j) );
            end
            A = (vn - vs)/hy + (ve - vw)/hx;
            D = nu/hy^2 * ( v(i,j+1)-2*v(i,j)+v(i,j-1) ) +...
                nu/hx^2 * ( v(i+1,j)-2*v(i,j)+v(i-1,j));
            
            vt(i,j)=v(i,j)+dt*(-A+D); 
        end
    end 
end

%-------------------CONCENTRATIONS BOUNDARY CONDITION----------------------
function [phi] = WallsBoundaryConditions(phi,nx,ny)

    phi(2:nx+1,1)=phi(2:nx+1,2);          
    phi(2:nx+1,ny+2)=phi(2:nx+1,ny+1);    
    phi(1,2:ny+1)=phi(2,2:ny+1);          
    phi(nx+2,2:ny+1)=phi(nx+1,2:ny+1); 

end

%----------------ADVECTION-DIFFUSION EQUATION FOR SPECIES------------------
function [phi]=TemporaryConcentrations(phi,u,v,dt,hx,hy,Diff,nx,ny)

    phio=phi;

    for i=2:nx+1
        for j=2:ny+1
            
            ue = u(i,j);    uw = u(i-1,j);
            vn = v(i,j);    vs = v(i,j-1);    
            
            phi(i,j)=phio(i,j)+dt*(-ue/(2*hx)*phio(i+1,j)+uw/(2*hx)*phio(i-1,j)+...
                                   -vn/(2*hy)*phio(i,j+1)+vs/(2*hy)*phio(i,j-1))+...
                               dt*Diff*((phio(i+1,j)-2*phio(i,j)+phio(i-1,j))/hx^2)+...
                               dt*Diff*((phio(i,j+1)-2*phio(i,j)+phio(i,j-1))/hy^2);                                              
        end
    end

end
%--------------------------REACTION STEP SOLVER----------------------------
function [phi] = BatchSolution(phi,nx,ny,dt_react,nu,r1,r2,r3)

    phio=phi;
    
    for i=2:nx+0
        for j = 2:ny+1

            phi(i,j) = phio(i,j)+dt_react*(nu(1)*r1(i,j)+nu(2)*r2(i,j)+nu(3)*r3(i,j));

        end
    end

end
