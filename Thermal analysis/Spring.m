R_E = 6348;%RADIUS OF EARTH
L = 1.5;%SIDE LENGTH OF A PLATE
A = L*L;%area of plates
DEN = 2700;%DENSITY OF ALUMINIUM
V = A*D;%volume of a plate
M = V*DEN;%mass of a plate
a = 0.3022;%alpha (radians)
DT = 60;%TIME STEP
T_I = 273;%INITIAL TEMPERATURE AT NODES

D = 0.01;%PLATE THICKNESS
C_AL = 896;%THERMAL CAPACITY OF ALUMINIUM
ABS = 0.2;%ABSORPTION COEFFICIENT
E = 0.05;%EMISSION COEFFICIENT 
K_AL = 229;%THERMAL CONDUCTIVITY
B = 5.66961e-8;%BOLTZMANN CONSTANT
E_0 = 1372;%SOLAR CONSTANT
E_IK = 0.2;%RADIATION VIEW FACTOR
BETA = 23.14;%SOLAR ALTITUDE IN SUMMER
R_G = 42164;%g=GEO RADIUS
phi(1) = 0;%initial rotation angle 
T1(1)=T_I;
T2(1)=T_I;
T3(1)=T_I;
T4(1)=T_I;
T5(1)=T_I;
T6(1)=T_I;
Q_C65 = zeros(1,5764);
 
for s = 1:1:5764%One step less to max
 
 %heat emission due to emission
 Q_E6(s)= E*B*2*A*(T6(s)^4);
 Q_E5(s)= E*B*2*A*(T5(s)^4);
 Q_E4(s)= E*B*2*A*(T4(s)^4);
 Q_E3(s)= E*B*2*A*(T3(s)^4);
 Q_E2(s)= E*B*2*A*(T2(s)^4);
 Q_E1(s)= E*B*2*A*(T1(s)^4);
  
 phi(s+1) = phi(s)+ (2*pi/1440);
 if phi(s+1) > 2*pi
    phi(s+1) = 0;
 else 
    phi(s+1) = phi(s+1);
 end
 
 %heat transfer due to solar radiation
 Q_S1(s) = E_0*ABS*(L^2)*cos(phi(s));
 if Q_S1(s)<0 || (phi(s)>(pi - a/2) && phi(s)<(pi + a/2))
     Q_S1(s) = 0;
 else 
     Q_S1(s) = Q_S1(s);
 end
Q_S2(s) = 0;
 Q_S3(s) = E_0*ABS*L^2*sin(phi(s));
 if Q_S3(s)<0 || (phi(s)>(pi - a/2) && phi(s)<(pi + a/2))
     Q_S3(s) = 0;
 else 
     Q_S3(s) = Q_S3(s);
 end
 
 Q_S4(s) = -(E_0*ABS*L^2*sin(phi(s)));
 if Q_S4(s)<0 || (phi(s)>(pi - a/2) && phi(s)<(pi + a/2))
     Q_S4(s) = 0;
 else 
     Q_S4(s) = Q_S4(s);
 end
 
 Q_S5(s) = 0;
 
 Q_S6(s) = -(E_0*ABS*L^2*cos(phi(s)));
 if Q_S6(s)<0 || (phi(s)>(pi - a/2) && phi(s)<(pi + a/2))
     Q_S6(s) = 0;
 else 
     Q_S6(s) = Q_S6(s);
 end
 
 %heat transfer due to conductance
 
 Q_C12(s) = K_AL*D*(T1(s)-T2(s));
 Q_C21(s) = -(Q_C12(s));
 Q_C13(s) = K_AL*D*(T1(s)-T3(s));
 Q_C31(s) = -(Q_C13(s));
 Q_C14(s) = K_AL*D*(T1(s)-T4(s));
 Q_C41(s) = -(Q_C14(s));
 Q_C15(s) = K_AL*D*(T1(s)-T5(s));
 Q_C51(s) = -(Q_C15(s));
 
 Q_C23(s) = K_AL*D*(T2(s)-T3(s));
 Q_C32(s) = -(Q_C23(s));
 Q_C24(s) = K_AL*D*(T2(s)-T4(s));
 Q_C42(s) = -(Q_C24(s));
 Q_C26(s) = K_AL*D*(T2(s)-T6(s));
 Q_C62(s) = -(Q_C26(s));
 
 Q_C46(s) = K_AL*D*(T4(s)-T6(s));
 Q_C64(s) = -(Q_C46(s));
 Q_C45(s) = K_AL*D*(T4(s)-T5(s));
 Q_C54(s) = -(Q_C45(s));
 Q_C56(s) = KAL*D*(T5(s)-T6(s));
 Q_C35(s) = K_AL*D*(T3(s)-T5(s));
 Q_C53(s) = -(Q_C35(s));
 Q_C36(s) = K_AL*D*(T3(s)-T6(s));
 Q_C63(s) = -(Q_C36(s));
 Q_C65(s) = -(Q_C56(s));
 


 %Heat absorption from other plates
 
 Q_A1(s) = (ABS*E_IK*E*B*A*(T2(s)^4 + T3(s)^4 + T4(s)^4 + T5(s)^4+ T6(s)^4));
 Q_A2(s) = (ABS*E_IK*E*B*A*(T1(s)^4 + T3(s)^4 + T4(s)^4 + T5(s)^4+ T6(s)^4));
 Q_A3(s) = (ABS*E_IK*E*B*A*(T2(s)^4 + T1(s)^4 + T4(s)^4 + T5(s)^4+ T6(s)^4));
 Q_A4(s) = (ABS*E_IK*E*B*A*(T2(s)^4 + T3(s)^4 + T1(s)^4 + T5(s)^4+ T6(s)^4));
 Q_A5(s) = (ABS*E_IK*E*B*A*(T2(s)^4 + T3(s)^4 + T4(s)^4 + T1(s)^4+ T6(s)^4));
 Q_A6(s) = (ABS*E_IK*E*B*A*(T2(s)^4 + T3(s)^4 + T4(s)^4 + T5(s)^4+ T1(s)^4));
 
 %Total Heat i.e. Absorbed-Emitted heat
 
 Q1(s) = Q_S1(s)+Q_A1(s)-Q_C12(s)-Q_C13(s)-Q_C14(s)-Q_C15(s)-Q_E1(s);
 Q2(s) = Q_S2(s)+Q_A2(s)-Q_C21(s)-Q_C23(s)-Q_C24(s)-Q_C26(s)-Q_E2(s);
 Q3(s) = Q_S3(s)+Q_A3(s)-Q_C31(s)-Q_C32(s)-Q_C36(s)-Q_C35(s)-Q_E3(s);
 Q4(s) = Q_S4(s)+Q_A4(s)-Q_C41(s)-Q_C42(s)-Q_C46(s)-Q_C45(s)-Q_E4(s);
 Q5(s) = Q_S5(s)+Q_A5(s)-Q_C51(s)-Q_C53(s)-Q_C54(s)-Q_C56(s)-Q_E5(s);
 Q6(s) = Q_S6(s)+Q_A6(s)-Q_C62(s)-Q_C63(s)-Q_C64(s)-Q_C65(s)-Q_E6(s);

 %temperature increment with change in phi per 60s
 
 T1(s+1) = (Q1(s)*60/(M*C_AL) + T1(s));
 T2(s+1) = (Q2(s)*60/(M*C_AL) + T2(s));
 T3(s+1) = (Q3(s)*60/(M*C_AL) + T3(s));
 T4(s+1) = (Q4(s)*60/(M*C_AL) + T4(s));
 T5(s+1) = (Q5(s)*60/(M*C_AL) + T5(s));
 T6(s+1) = (Q6(s)*60/(M*C_AL) + T6(s));
end
figure(1);
legend('T6','T5','T4','T3','T2','T1');
plot(T6,'r')
hold on
plot(T5,'g')
plot(T4,'m')
plot(T3,'y')
plot(T2,'b')
plot(T1,'k') 
