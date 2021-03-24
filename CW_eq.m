function tmp_eq = CW_eq(t,X_0)
%CW_eq This function encompasses the CW equations
% The below equations are created assuming 2 assumptions
% 1.If (x,y,z) are much smaller than r_c, where r_c is radius of cheif satelleite's orbit
% 2.The Chief satellite's orbit is circular
%
%   ddot(x) - 2ndot(y) - 3n^2x = 0
%   ddot(y) + 2ndot(x) = 0
%   ddot(z) + n^2z = 0

%% Creating the symbolic variables and functions
syms x(t) y(t) z(t) %Variables that depend on time
syms n %Variables that do not depend on time

%Below the CW equations are defined.

eq1 = diff(x,2) - 2*n*diff(y,1) - 3*n^2*x == 0;
eq2 = diff(y,2) + 2*n*diff(x,1) == 0;
eq3 = diff(z,2) + n^2*z == 0;

[eq_sym,var_sym] = odeToVectorField(eq1,eq2,eq3);

tmp_eq = matlabFunction(eq_sym,'vars',{'t','r'});









end

