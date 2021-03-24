% Defining the system
U_s = 300; % 
m = 500; % Defining the mass of the body
K_M = 10; %
K_R = 0.25 ; %
v_s = sqrt(K_M / K_R * U_s ) ; %
u = U_s * ones (200 ,1) ; %
v = zeros (200 ,1) ; %
a = zeros (200 ,1) ; %
x = zeros (200 ,1) ; %
% Defining velocity and acceleration
v ( 1 ) = 10; %
a ( 1 ) = 6; %
dt = 1; %
% Calculations
for t = 2:199  % Taking time in the range of 2s to 199s
v ( t ) = v (t-1) + K_R/m * v (t -1)^2 + K_M/m*u (t)*dt ; % calculating velocity as a function of time
a ( t ) = ( v ( t ) - v ( t - 1) / dt ) ; % calculating acceleration as a function of time
end % We end the for loop here
% Graphs
subplot (311) %
t =1:200; %
plot ( t , u ) ; % Plot between Time and initial velocity
subplot (312) %
t =1:200; %
plot ( t , v ) ; % Plot between time and final velocity
subplot (313) %
t =1:200; %
plot ( t , a ) ; % Plot between time and acceleration