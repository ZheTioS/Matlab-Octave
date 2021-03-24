s = tf ( ’ s ’ ) ;
I = 64;
G_R = 10; % P-Controller
G_R = s ; % I-Controller
G_R = 10/s; % D-Controller
G_R = 10+s; % PD-Controller
G_R = 10+10/s+s; % PID-Controller
G_R = 10+10/s; % PI-Controller
G_R = 10/s+s; % ID-Controller
G_S = ( 1 / I ) / s ^2;
G_0 = G_S*G_R; % Open Loop transfer function
G = G_0/1+G_0; % Closed loop transfer function
figure ( 1 )
impulse (G)
figure ( 2 )
step (G)
