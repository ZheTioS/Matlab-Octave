x_0 = 0; %m
y_0 = 100; %m
z_0 = 0; %m

a = 7000; %km
mu = 398600.4418; %km^3/s

v_xin = 0.05785; %m/s
v_yin = 0; %m/s
v_zin = 0.11570; %m/s
r = [x_0; y_0; z_0];
v = [v_xin; v_yin; v_zin];
[a,eMag,i,O,o,nu,truLon,argLat,lonPer,p] = rv2orb(r,v,mu);
i = rad2deg(i)
O = rad2deg(O)
o = rad2deg(o)