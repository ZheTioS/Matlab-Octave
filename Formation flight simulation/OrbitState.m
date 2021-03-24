function [State] = OrbitState(t,S,u,K,Body_Radius,J_num,flag)
    % find new position magnitude at each time step
    r_new = norm(S(1:3));
    % verify to perturb the orbit or not
    switch flag
        case 1
            % Initialize symbolic variables
            syms u_param r r_e phi k a x y z;
            A_d_old = 0; % initialize old value for disturbance acceleration
            % obtain A_d expressions 1:K
            for i=1:K
                [A_d_Expression,J] = Get_Disturbance(i,u_param,r,r_e,phi,k,a,x,y,z); % call the Get_Disturbance() function
                A_d_Expression = A_d_Expression + A_d_old; % add the previous value to the new value
                A_d_old = A_d_Expression; % previous value, becomes new value
            end            
            % substitute constant numerical values
            A_d_Expression = subs(A_d_Expression,r_e,Body_Radius); % sub for equatorial radius
            A_d_Expression = subs(A_d_Expression,u_param,u); % sub for gravitational parameter
            % substitute numerical values which are not constant into the expression
            A_d_Expression = subs(A_d_Expression,x,S(1)); % sub for x
            A_d_Expression = subs(A_d_Expression,y,S(2)); % sub for y
            A_d_Expression = subs(A_d_Expression,z,S(3)); % sub for z
            A_d_Expression = subs(A_d_Expression,r,r_new); % sub for r magnitude
            for i=1:K
                A_d_Expression = subs(A_d_Expression,J(1,i),J_num(i)); % how i do this?
            end
            A_d = double(vpa(A_d_Expression,6)); % simplify numerical values out to 6 sigfigs
            
        otherwise
            A_d = [0.0;0.0;0.0]; % set disturbance acceleration to zero
    end
        % Output State from orbital state-space form
    State = [S(4)
             S(5)
             S(6)
             (-u/r_new^3)*S(1) + A_d(1)
             (-u/r_new^3)*S(2) + A_d(2)
             (-u/r_new^3)*S(3) + A_d(3)];
end