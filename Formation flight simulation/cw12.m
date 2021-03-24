clear all
clc
%% Question 1.2
%% Defining the given information
x_0 = 0; %m
y_0 = 0; %m
z_0 = 0; %m

a = 7000; %km
nu = 398600.4418; %km^3/s

v_xin = 0; %m/s
v_yin = 0.1; %m/s
v_zin = 0; %m/s

%% Defining symbols for the CW equation.

syms x(t) y(t) z(t) %Variables that depend on time
syms n %Variables that do not depend on time
n = nu/a^3; % defining the value of n
%Below the CW equations are defined.

eq1 = diff(x,2) - 2*n*diff(y,1) - 3*n^2*x == 0;
eq2 = diff(y,2) + 2*n*diff(x,1) == 0;
eq3 = diff(z,2) + n^2*z == 0;

[sym_eq,sym_var] = odeToVectorField(eq1,eq2,eq3);

tmp_eq = matlabFunction(sym_eq,'vars',{'t','Y'});
%% Now using Runge-kutta method

X_0 = [x_0 v_xin y_0 v_yin z_0 v_zin]; %storing initial values into 1 variable
tspan= 0:6000; % Running RK for a time steps upto 100
[T , pos] = ode45(tmp_eq,tspan,X_0); %using ODE45 which uses RK4 method.

x_2 = pos(:,1);
vx_2 = pos(:,2);
y_2 = pos(:,3);
vy_2 = pos(:,4);
z_2 = pos(:,5);
vz_2 = pos(:,6);

figure(2)

plot3(x_2,y_2,z_2)