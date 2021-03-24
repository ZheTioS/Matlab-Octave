syms x(t) y(t) z(t) %Variables that depend on time
syms n %Variables that do not depend on time
nu = 3986000;
a = 7000;
n = nu/a^3;
%Below the CW equations are defined.

eq1 = diff(x,2) - 2*n*diff(y,1) - 3*n^2*x == 0;
eq2 = diff(y,2) + 2*n*diff(x,1) == 0;
eq3 = diff(z,2) + n^2*z == 0;

[eq_sym,var_sym] = odeToVectorField(eq1,eq2,eq3)

tmp_eq = matlabFunction(eq_sym,'vars',{'t','Y'})

init = [0 0.5 10 0 0 0.5]
tspan=[0:20];
sol = ode45(tmp_eq,tspan,init)