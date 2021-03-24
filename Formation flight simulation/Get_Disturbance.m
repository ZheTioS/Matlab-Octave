function [A_d_Expression,J] = Get_Disturbance(K,u_param,r,r_e,phi,k,a,x,y,z)
% INPUT = all symbolic variables
% OUTPUT = Symbolic expression in Cartesian coordinates for disturbance acceleration, symbolic J terms
    
J = sym('J',[1,K]); % create J according to polynomial expansion coefficient, J1,J2, or J3,... 
    % Setup expression for Legendre Polynomials
    f = (a^2 - 1)^K;
    PK = (1/(2^(K)*factorial(K)))*diff(f,a,K);
    % substitute eigenfunction
    a = cos(phi);   
    % substitute for PK
    PK = subs(PK);
    % g expression from infinite series   
    gK = J(K)*(r_e/r)^K*PK;
    gK = subs(gK);
    UK = -(u_param/r)*symsum(gK,k,K,K);
     % Take the negative gradient
    A_d = -[diff(UK,r);diff(UK,phi)];
    phi = acos(z/(sqrt(x^2 + y^2 + z^2)));
    A_d = subs(A_d);
    dPhi = [diff(phi,x);diff(phi,y);diff(phi,z)];
    % cartesian conversion
    A_d_Expression = [A_d(1)*x/r + A_d(2)*dPhi(1)
                  A_d(1)*y/r + A_d(2)*dPhi(2)
                  A_d(1)*z/r + A_d(2)*dPhi(3)];
    A_d_Expression = subs(A_d_Expression,x^2 + y^2 + z^2,r^2);
end