s = tf('s');
G_s = (1)/(s^2);
PI = 1+ 1/s;
PD = 1+s;
G = G_s*PD;
rlocus(G)