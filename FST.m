%% FST-Formation Flight

clear all; 
clc;
%% Question 1.1
%% Defining the given information
x0 = 0 %m
y0 = 100 %m
z0 = 0 %m

a = 7000 %km
nu = 398600.4418 %km^3/s

v_xi = 0.05785 %m/s
v_yi = 0 %m/s
v_zi = 0.11570 %m/s

%% The CW equation

% The below equations are created assuming 2 assumptions

% 1. If (x,y,z) are much smaller than r_c, where r_c is radius of cheif
% satelleite's orbit.

% 2. The Chief satellite's orbit is circular.
%
%   ddot(x) - 2ndot(y) - 3n^2x = 0
%   ddot(y) + 2ndot(x) = 0
%   ddot(z) + n^2z = 0
%% Defining symbols for the CW equation.
syms x(t) y(t) z(t) %Variables that depend on time
syms n %Variables that do not depend on time
n = sqrt(nu/a^3) % defining the value of n
%Below the CW equations are defined.

eq1 = diff(x,2) - 2*n*diff(y,1) - 3*n^2*diff(x,0) == 0
eq2 = diff(y,2) + 2*n*diff(x,1) == 0
eq3 = diff(z,2) + n^2*diff(z,0) == 0

[sym_eq,sym_var] = odeToVectorField(eq1,eq2,eq3)

tmp_eq = matlabFunction(sym_eq,'vars',{'t','Y'})
%% Now using Runge-kutta method
X_0 = [x0 v_xi y0 v_yi z0 v_zi] %storing initial values into 1 variable
tspan= 0:6000; % Running RK for a time steps upto 6000
[~ , Pos] = ode45(tmp_eq,tspan,X_0) %using ODE45 which uses RK4 method.

x_1 = Pos(:,1)
vx_1 = Pos(:,2)
y_1 = Pos(:,3)
vy_1 = Pos(:,4)
z_1 = Pos(:,5)
vz_1 = Pos(:,6)

figure(1)
plot3(x_1,y_1,z_1)
grid on
xlabel('X(km)')
ylabel('Y(km)')
zlabel('Z(km)')

%% Question 1.2
%% Defining the given information
x_0 = 0 %m
y_0 = 0 %m
z_0 = 0 %m

a = 7000 %km
nu = 398600.4418 %km^3/s

v_xin = 0; %m/s
v_yin = 0.1; %m/s
v_zin = 0; %m/s

%% Defining symbols for the CW equation.

syms x(t) y(t) z(t) %Variables that depend on time
syms n %Variables that do not depend on time
n = sqrt(nu/a^3) % defining the value of n
%Below the CW equations are defined.

eq1 = diff(x,2) - 2*n*diff(y,1) - 3*n^2*x == 0
eq2 = diff(y,2) + 2*n*diff(x,1) == 0
eq3 = diff(z,2) + n^2*z == 0

[sym_eq,sym_var] = odeToVectorField(eq1,eq2,eq3)

tmp_eq = matlabFunction(sym_eq,'vars',{'t','Y'})
%% Now using Runge-kutta method

X_0 = [x_0 v_xin y_0 v_yin z_0 v_zin]; %storing initial values into 1 variable
tspan= 0:60000; % Running RK for a time steps
[~, pos] = ode45(tmp_eq,tspan,X_0); %using ODE45 which uses RK4 method.

x_2 = pos(:,1);
vx_2 = pos(:,2);
y_2 = pos(:,3);
vy_2 = pos(:,4);
z_2 = pos(:,5);
vz_2 = pos(:,6);

figure(2)

plot3(x_2,y_2,z_2)
grid on
xlabel('X(km)');
ylabel('Y(km)');
zlabel('Z(km)');

% It is observed that the deputy moves away from the Chief
%% Question 2
%Given Orbital elements for chief are
% M_0 = deg2rad(0);
% e = deg2rad(0);
% i = deg2rad(98);
% RAAN = deg2rad(0);
% arg_peri = deg2rad(45);
% a = 7000;
% 
% %% Conversion to cartesian co-ordinates [1]
% [r, v] = orb2rv(a,e,i,RAAN,arg_peri,M_0,0,0,0,nu);

%% Setting Initial values
x_c = 4.949747468305833 * 10^3
y_c =  -0.688871704133356 * 10^3
z_c = 4.901576866197694 * 10^3 
vx_c =  -5.335862495551078 
vy_c = -0.742608529802357  
vz_c = 5.283934248539943 

%% Equations to define orbit of Chief

Position = [x_c;y_c;z_c]
Velocity = [vx_c;vy_c;vz_c]
%Defining Earth as per [2]
        u = 3.986004415e5; % km^3/s^2
        Body = 'Earth'
        map = imread('Images/Earth_map.jpg'); % read planet map
        map = imrotate(map,180); % rotate upright
        Body_Radius = 6371; % km
        Body_Mass = 5.97219e24; % kg
        Dist_to_Sun = 149.6e6; % km
        J_num = [0;0.0010826269;-0.0000025323;-0.0000016204]; % zonal harmonic coefficients
% find angular momentum km^2/s
hVec = cross(Position,Velocity); % vector
h = norm(hVec); % magnitude
e0 = cross(Position,hVec)
em0=norm(e0);
% find eccentricity
r = norm(Position); % magnitude
eVec = (cross(Velocity,hVec)/u) - (Position/r); % vector
e = norm(eVec); % magnitude
% find semi-major axis in km
a = h^2/(u*(1-e^2))

% Place initial system in matrix Form
sys_0 = [Position(1);Position(2);Position(3);Velocity(1);Velocity(2);Velocity(3)];

Norbits = 3 % number of orbits to propagate (only matters in a perturbed circular or elliptical orbit)

% Create a function to handle the calculation of the time of flight
[T,Orbit] = Get_TOF(a,u,e,Body_Mass,Dist_to_Sun,Norbits);
% develop a time interval
dt = 50; % time step
tspan = 0:dt:T; % this is in seconds
Sim = 'Unperturbed'; % specify how to run the sim ('Perturbed', 'Unperturbed')
% Simulation using ode45 and OrbitState function
options = odeset('RelTol',1e-8,'AbsTol',1e-9); % specify ODE45 tolerances
switch Sim
    case 'Unperturbed'
        flag = 0 % Chosen to unperturb orbit
        K = []; % K is empty
        [t,S] = ode45(@(t,S) OrbitState(t,S,u,K,Body_Radius,J_num,flag), tspan, sys_0, options); % propagate
        Unperturbed_States = S;
    case 'Perturbed'
        flag = 1; % Chosen to perturb orbit
        K = numel(J_num) % obtain the number of elements in J, this represents the level of fidelity available to find disturbance acceleration
        [t,S] = ode45(@(t,S) OrbitState(t,S,u,K,Body_Radius,J_num,flag), tspan, sys_0, options); % propagate
        Perturbed_States = S;
end
% Separate my states into x,y,z components for chief
States_X_C = S(:,1); % all x
States_Y_C = S(:,2); % all y
States_Z_C = S(:,3); % all z
States_Xdot_C = S(:,4); % all xdot
States_Ydot_C = S(:,5); % all ydot
States_Zdot_C = S(:,6); % all zdot

%% 2.1 For deputy

%Adding the relative values from 1.1
Position = [x_0+x_c;y_0+y_c;z_0+z_c]
Velocity = [v_xi+vx_c;v_yi+vy_c;v_zi+vz_c]
%Defining Earth as per [2]
        u = 3.986004415e5; % km^3/s^2
        Body = 'Earth';
        map = imread('Images/Earth_map.jpg'); % read planet map
        map = imrotate(map,180); % rotate upright
        Body_Radius = 6371; % km
        Body_Mass = 5.97219e24; % kg
        Dist_to_Sun = 149.6e6; % km
        J_num = [0;0.0010826269;-0.0000025323;-0.0000016204]; % zonal harmonic coefficients
% find angular momentum km^2/s
hVec= cross(Position,Velocity); % vector
h= norm(hVec); % magnitude
% find eccentricity
r = norm(Position); % magnitude
eVec = (cross(Velocity,hVec)/u) - (Position/r) % vector
e = norm(eVec); % magnitude
% find semi-major axis in km
a = h^2/(u*(1-e^2))

% Place initial system in matrix Form
sys_0 = [Position(1);Position(2);Position(3);Velocity(1);Velocity(2);Velocity(3)];

Norbits = 3 % number of orbits to propagate (only matters in a perturbed circular or elliptical orbit)

% Create a function to handle the calculation of the time of flight
[T,Orbit] = Get_TOF(a,u,e,Body_Mass,Dist_to_Sun,Norbits);
% develop a time interval
dt = 50; % time step
tspan = 0:dt:T; % this is in seconds
Sim = 'Unperturbed'; % specify how to run the sim ('Perturbed', 'Unperturbed')
% Simulation using ode45 and OrbitState function
options = odeset('RelTol',1e-8,'AbsTol',1e-9); % specify ODE45 tolerances
switch Sim
    case 'Unperturbed'
        flag = 0; % Chosen to unperturb orbit
        K = []; % K is empty
        [t,S] = ode45(@(t,S) OrbitState(t,S,u,K,Body_Radius,J_num,flag), tspan, sys_0, options); % propagate
        Unperturbed_States = S;
    case 'Perturbed'
        flag = 1; % Chosen to perturb orbit
        K = numel(J_num) % obtain the number of elements in J, this represents the level of fidelity available to find disturbance acceleration
        [t,S] = ode45(@(t,S) OrbitState(t,S,u,K,Body_Radius,J_num,flag), tspan, sys_0, options); % propagate
        Perturbed_States = S;
end
% Separate my states into x,y,z components for deputy
States_X_D = S(:,1) % all x
States_Y_D = S(:,2) % all y
States_Z_D = S(:,3) % all z
States_Xdot_D = S(:,4) % all xdot
States_Ydot_D = S(:,5) % all ydot
States_Zdot_D = S(:,6) % all zdot

%% Develop the figure plot for Chief and Deputy Trajectory
figure(3);
plot3(States_X_C,States_Y_C,States_Z_C,'r');
grid on;
title_form = sprintf('Chief and Deputy(1.1) Trajectory'); % make the title unique to the orbit body and conic section
title(title_form);
xlabel('X(km)');
ylabel('Y(km)');
zlabel('Z(km)');
 [x,y,z] = sphere; % create a sphere
 Cord = [0 0 0 Body_Radius]; % use the equatorial radius of orbit body as a coordinate paramater
hold on;
plot3(States_X_D,States_Y_D,States_Z_D,'b');
 Obj = surf(x*Cord(1,4),y*Cord(1,4),z*Cord(1,4)); % create a surface plot unique to the equatorial radius
 set(Obj,'edgecolor','none');
 daspect([1 1 1]); % set the aspect ratio
 warp(x*Cord(1,4),y*Cord(1,4),z*Cord(1,4),map);
hold off;
legend('Chief','Deputy')

%% 2.2 Deputy
% Using the Deputy case from 1.2
Position = [x_0+x_c;y_0+y_c;z_0+z_c]; 
Velocity = [v_xin+vx_c;v_yin+vy_c;v_zin+vz_c];
%Defining Earth as per [2]
        u = 3.986004415e5; % km^3/s^2
        Body = 'Earth';
        map = imread('Images/Earth_map.jpg'); % read planet map
        map = imrotate(map,180); % rotate upright
        Body_Radius = 6371; % km
        Body_Mass = 5.97219e24; % kg
        Dist_to_Sun = 149.6e6; % km
        J_num = [0;0.0010826269;-0.0000025323;-0.0000016204]; % zonal harmonic coefficients
% find angular momentum km^2/s
hVec = cross(Position,Velocity); % vector
h = norm(hVec); % magnitude
% find eccentricity
r = norm(Position); % magnitude
eVec = (cross(Velocity,hVec)/u) - (Position/r); % vector
e = norm(eVec); % magnitude
% find semi-major axis in km
a = h^2/(u*(1-e^2));

% Place initial system in matrix Form
sys_0 = [Position(1);Position(2);Position(3);Velocity(1);Velocity(2);Velocity(3)];

Norbits = 3; % number of orbits to propagate (only matters in a perturbed circular or elliptical orbit)

% Create a function to handle the calculation of the time of flight
[T,Orbit] = Get_TOF(a,u,e,Body_Mass,Dist_to_Sun,Norbits);
% develop a time interval
dt = 50; % time step
tspan = 0:dt:T; % this is in seconds
Sim = 'Unperturbed'; % specify how to run the sim ('Perturbed', 'Unperturbed')
% Simulation using ode45 and OrbitState function
options = odeset('RelTol',1e-8,'AbsTol',1e-9); % specify ODE45 tolerances
switch Sim
    case 'Unperturbed'
        flag = 0; % Chosen to unperturb orbit
        K = []; % K is empty
        [t,S] = ode45(@(t,S) OrbitState(t,S,u,K,Body_Radius,J_num,flag), tspan, sys_0, options); % propagate
        Unperturbed_States = S;
    case 'Perturbed'
        flag = 1; % Chosen to perturb orbit
        K = numel(J_num) % obtain the number of elements in J, this represents the level of fidelity available to find disturbance acceleration
        [t,S] = ode45(@(t,S) OrbitState(t,S,u,K,Body_Radius,J_num,flag), tspan, sys_0, options); % propagate
        Perturbed_States = S;
end
% Separate my states into x,y,z components for deputy
States_X_D = S(:,1) % all x
States_Y_D = S(:,2) % all y
States_Z_D = S(:,3) % all z
States_Xdot_D = S(:,4) % all xdot
States_Ydot_D = S(:,5) % all ydot
States_Zdot_D = S(:,6) % all zdot

%% Develop the figure plot for Chief and Deputy Trajectory

figure(4);
plot3(States_X_C,States_Y_C,States_Z_C,'r');
grid on;
title_form = sprintf('Chief and Deputy(1.2) Trajectory');
title(title_form);
xlabel('X(km)');
ylabel('Y(km)');
zlabel('Z(km)');
 [x,y,z] = sphere; % create a sphere
 Cord = [0 0 0 Body_Radius]; % use the equatorial radius of orbit body as a coordinate paramater
hold on;
plot3(States_X_D,States_Y_D,States_Z_D,'b');
 Obj = surf(x*Cord(1,4),y*Cord(1,4),z*Cord(1,4)); % create a surface plot unique to the equatorial radius
 set(Obj,'edgecolor','none');
 daspect([1 1 1]); % set the aspect ratio
 warp(x*Cord(1,4),y*Cord(1,4),z*Cord(1,4),map);
hold off;
legend('Chief','Deputy')

%% Question 3

Position = [-x0+x_c;-y0+y_c;-z0+z_c]; 
Velocity = [-v_xi+vx_c;-v_yi+vy_c;-v_zi+vz_c];
%Defining Earth as per [2]
        u = 3.986004415e5; % km^3/s^2
        Body = 'Earth';
        map = imread('Images/Earth_map.jpg'); % read planet map
        map = imrotate(map,180); % rotate upright
        Body_Radius = 6371; % km
        Body_Mass = 5.97219e24; % kg
        Dist_to_Sun = 149.6e6; % km
        J_num = [0;0.0010826269;-0.0000025323;-0.0000016204]; % zonal harmonic coefficients
% find angular momentum km^2/s
hVec = cross(Position,Velocity); % vector
h = norm(hVec); % magnitude
% find eccentricity
r = norm(Position); % magnitude
eVec = (cross(Velocity,hVec)/u) - (Position/r); % vector
e = norm(eVec); % magnitude
% find semi-major axis in km
a = h^2/(u*(1-e^2));

% Place initial system in matrix Form
sys_0 = [Position(1);Position(2);Position(3);Velocity(1);Velocity(2);Velocity(3)];

Norbits = 3; % number of orbits to propagate (only matters in a perturbed circular or elliptical orbit)

% Create a function to handle the calculation of the time of flight
[T,Orbit] = Get_TOF(a,u,e,Body_Mass,Dist_to_Sun,Norbits);
% develop a time interval
dt = 50; % time step
tspan = 0:dt:T; % this is in seconds
Sim = 'Unperturbed'; % specify how to run the sim ('Perturbed', 'Unperturbed')
% Simulation using ode45 and OrbitState function
options = odeset('RelTol',1e-8,'AbsTol',1e-9); % specify ODE45 tolerances
switch Sim
    case 'Unperturbed'
        flag = 0; % Chosen to unperturb orbit
        K = []; % K is empty
        [t,S] = ode45(@(t,S) OrbitState(t,S,u,K,Body_Radius,J_num,flag), tspan, sys_0, options); % propagate
        Unperturbed_States = S;
    case 'Perturbed'
        flag = 1; % Chosen to perturb orbit
        K = numel(J_num) % obtain the number of elements in J, this represents the level of fidelity available to find disturbance acceleration
        [t,S] = ode45(@(t,S) OrbitState(t,S,u,K,Body_Radius,J_num,flag), tspan, sys_0, options); % propagate
        Perturbed_States = S;
end
% Separate my states into x,y,z components for deputy
States_X_D = S(:,1) % all x
States_Y_D = S(:,2) % all y
States_Z_D = S(:,3) % all z
States_Xdot_D = S(:,4) % all xdot
States_Ydot_D = S(:,5) % all ydot
States_Zdot_D = S(:,6) % all zdot

%% Develop the figure plot for Chief and Deputy Trajectory

figure;
grid on;
title_form = sprintf('Chief and Deputy Relative Trajectory'); 
title(title_form);
xlabel('X(km)');
ylabel('Y(km)');
zlabel('Z(km)');
hold on;
plot3(States_X_D,States_Y_D,States_Z_D,'y');
plot3(x_1+Body_Radius,y_1+Body_Radius,z_1+Body_Radius,'r');
hold off;
legend('Relative','1.1')


%% Refrences
%1.Darin Koblick (2020). Convert Keplerian Orbital Elements to a State Vector 
%  (https://www.mathworks.com/matlabcentral/fileexchange/35455-convert-keplerian-orbital-
%  elements-to-a-state-vector), MATLAB Central File Exchange. Retrieved January 27, 2020.
%2.Richard DeMark (2020). Orbit Trajectory Propagation 
%  (https://www.mathworks.com/matlabcentral/fileexchange/69439-orbit-trajectory-propagation),
%  MATLAB Central File Exchange. Retrieved January 29, 2020.