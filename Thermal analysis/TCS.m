clc
%Program to perform Temperature Analysis
%% Given Values are
Cth = 896;% Thermal Conductivity for Al
RHO = 2.7*10^(3); % Density of Al
ABS = 0.3; %Absorbtion Co-efficient for Al
EAL = 0.1; %Emission Co-efficient of Al
KAL = 229; %Thermal Conductivity of Al
SIG = 5.9961*10^(-9); %Boltzmann Co-efficient
E0 = 1372; %Solar Constant
e_ik = 0.2; %Radiation view factor
dt = 60; %Time step
NOS = 5; %Number of Orbits for simulation
T_in = 290; %Inital temperature for Nodes
L= 2; %Side Length
d = 0.02; %Plate Thickeness
BETA = 23.14;%Solar attitude in summer
rGEO= 42164;%Radius GEO taking earths radius into account
rE = 6378; %Radius of Earth
T1=zeros(1,7200);
T2=zeros(1,7200);
T3=zeros(1,7200);
T4=zeros(1,7200);
T5=zeros(1,7200);
T6=zeros(1,7200);
T1(1) = T_in;
T2(1) = T_in;
T3(1) = T_in;
T4(1) = T_in;
T5(1) = T_in;
T6(1) = T_in;
D = L*d;
m = RHO*L^2*d;
   Qem1 = zeros(1,7200);
   Qem2 = zeros(1,7200);
   Qem3 = zeros(1,7200);
   Qem4 = zeros(1,7200);
   Qem5 = zeros(1,7200);
   Qem6 = zeros(1,7200);
%% Projected Area of Nodes
%% Spring

A_spr = 2*asind(rE/rGEO); 
A1_sp=zeros(1,7200);
A2_sp=zeros(1,7200);
A3_sp=zeros(1,7200);
A4_sp=zeros(1,7200);
A5_sp=zeros(1,7200);
A6_sp=zeros(1,7200);

for PSI= 0:7
    i = 1+PSI;
    if((PSI >(180-(A_spr/2))) && (PSI <(180+(A_spr/2))))
        A1_sp(i) = 0;
        A3_sp(i) = 0;
        A4_sp(i) = 0;
        A6_sp(i) = 0; 
    else
        A1_sp(i) = L^2*cosd(PSI);
        A3_sp(i) = L^2*sind(PSI);
        A4_sp(i) = -(L^2*sind(PSI));
        A6_sp(i) = -(L^2*cosd(PSI));
    end
   
end
%% Summer

A1_sm=zeros(1,7200);
A2_sm=zeros(1,7200);
A3_sm=zeros(1,7200);
A4_sm=zeros(1,7200);
A5_sm=zeros(1,7200);
A6_sm=zeros(1,7200);
for PSI= 0:7200
    i = 1+PSI;
        A1_sm(i) = L^2*cosd(PSI)*sind(BETA);
        A2_sm(i) = L^2 * sind(BETA);
        A3_sm(i) = L^2*sind(PSI)*cosd(BETA);
        A4_sm(i) = -(L^2*sind(PSI)*cosd(BETA));
        A6_sm(i) = -(L^2*cosd(PSI)*cosd(BETA));
 
end

%% Calculations
h = 0.7; 
n = 7200; % No.of minutes in 5 days.
%% Spring

F = @(x,y)((x*(dt/(m*Cth)))+y);

for i = 1:n
%% Solar influence
QS_1 = zeros(1,7200);
QS_2 = zeros(1,7200);
QS_3 = zeros(1,7200);
QS_4 = zeros(1,7200);
QS_5 = zeros(1,7200);
QS_6 = zeros(1,7200);

for q = 1:7200
   if(A1_sp(q)~=0)
       QS_1(q)= E0*ABS*A1_sp(q);
   end
   if(A2_sp(q)~=0)
       QS_2(q)=E0*ABS*A2_sp(q);
   end
   if(A3_sp(q)~=0)
       QS_3(q)= E0*ABS*A3_sp(q);
   end
   if(A4_sp(q)~=0)
       QS_4(q)= E0*ABS*A4_sp(q);
   end
   if(A6_sp(q)>=0)
       QS_6(q)= E0*ABS*A6_sp(q);
   end
   
end
% plot(QS_1) 
% hold on 
% plot(QS_3)

%% Emission
   %Conductivity due to emission
   
   Qem1(i) = EAL*SIG*2*L^2*T1(i)^4; 
   Qem2(i) = EAL*SIG*2*L^2*T2(i)^4;
   Qem3(i) = EAL*SIG*2*L^2*T3(i)^4;
   Qem4(i) = EAL*SIG*2*L^2*T4(i)^4;
   Qem5(i) = EAL*SIG*2*L^2*T5(i)^4;
   Qem6(i) = EAL*SIG*2*L^2*T6(i)^4;
  
%% Thermal Conductivity b/w plates
 %Face 1
 
 QTC_12(i) = KAL*D*(T1(i)-T2(i))/L;
 QTC_13(i) = KAL*D*(T1(i)-T3(i))/L;
 QTC_14(i) = KAL*D*(T1(i)-T4(i))/L;
 QTC_15(i) = KAL*D*(T1(i)-T5(i))/L;
 
 %Face 2
 
 QTC_21(i) = -(QTC_12(i));
 QTC_23(i) = KAL*D*(T2(i)-T3(i))/L;
 QTC_24(i) = KAL*D*(T2(i)-T4(i))/L;
 QTC_26(i) = KAL*D*(T2(i)-T6(i))/L;
 
 %Face 3
 
 QTC_31(i) = -(QTC_13(i));
 QTC_32(i) = -(QTC_23(i));
 QTC_35(i) = KAL*D*(T3(i)-T5(i))/L;
 QTC_36(i) = KAL*D*(T3(i)-T6(i))/L;
  
 %Face 4
 
 QTC_41(i) = -(QTC_14(i));
 QTC_42(i) = -(QTC_24(i));
 QTC_46(i) = KAL*D*(T4(i)-T6(i))/L;
 QTC_45(i) = KAL*D*(T4(i)-T5(i))/L;
 
 %Face 5
 
 QTC_51(i) = -(QTC_15(i));
 QTC_54(i) = -(QTC_45(i));
 QTC_53(i) = -(QTC_35(i));
 QTC_56(i) = KAL*D*(T5(i)-T6(i))/L;
 
 %Face 6
 
 QTC_62(i) = -(QTC_26(i));
 QTC_64(i) = -(QTC_46(i));
 QTC_63(i) = -(QTC_36(i));
 QTC_65(i) = -(QTC_56(i));

 %% Radiation from other plates
 Q_A1(i) = (ABS*e_ik*EAL*SIG*L^2*(T2(i)^4 + T3(i)^4 + T4(i)^4 + T5(i)^4+ T6(i)^4));
 Q_A2(i) = (ABS*e_ik*EAL*SIG*L^2*(T1(i)^4 + T3(i)^4 + T4(i)^4 + T5(i)^4+ T6(i)^4));
 Q_A3(i) = (ABS*e_ik*EAL*SIG*L^2*(T2(i)^4 + T1(i)^4 + T4(i)^4 + T5(i)^4+ T6(i)^4));
 Q_A4(i) = (ABS*e_ik*EAL*SIG*L^2*(T2(i)^4 + T3(i)^4 + T1(i)^4 + T5(i)^4+ T6(i)^4));
 Q_A5(i) = (ABS*e_ik*EAL*SIG*L^2*(T2(i)^4 + T3(i)^4 + T4(i)^4 + T1(i)^4+ T6(i)^4));
 Q_A6(i) = (ABS*e_ik*EAL*SIG*L^2*(T2(i)^4 + T3(i)^4 + T4(i)^4 + T5(i)^4+ T1(i)^4));

 %% Calculating the overall conductivity on the system
 
 Q1(i) = -QS_1(i)+Q_A1(i)-QTC_12(i)-QTC_13(i)-QTC_14(i)-QTC_15(i)-Qem1(i);
 Q2(i) = -QS_2(i)+Q_A2(i)-QTC_21(i)-QTC_23(i)-QTC_24(i)-QTC_26(i)-Qem2(i);
 Q3(i) = -QS_3(i)+Q_A3(i)-QTC_31(i)-QTC_32(i)-QTC_36(i)-QTC_35(i)-Qem3(i);
 Q4(i) = -QS_4(i)+Q_A4(i)-QTC_41(i)-QTC_42(i)-QTC_46(i)-QTC_45(i)-Qem4(i);
 Q5(i) = -QS_5(i)+Q_A5(i)-QTC_51(i)-QTC_53(i)-QTC_54(i)-QTC_56(i)-Qem5(i);
 Q6(i) = -QS_6(i)+Q_A6(i)-QTC_62(i)-QTC_63(i)-QTC_64(i)-QTC_65(i)-Qem6(i);
 
 %% Calculation of temperature
     % Apply Runge Kutta Formulas to find next value of T
        %For Face 1
        k1_1 = h*F(Q1(i), T1(i)); 
        k2_1 = h*F(Q1(i) + 0.5*h, T1(i) + 0.5*k1_1); 
        k3_1 = h*F(Q1(i) + 0.5*h, T1(i) + 0.5*k2_1); 
        k4_1 = h*F(Q1(i) + h, T1(i) + k3_1);
        T1(i+1) = T1(i) + (1/6)*(k1_1 + 2*k2_1 + 2*k3_1 + k4_1);
        %for Face 2
        k1_2 = h*F(Q2(i), T2(i)); 
        k2_2 = h*F(Q2(i) + 0.5*h, T2(i) + 0.5*k1_2); 
        k3_2 = h*F(Q2(i) + 0.5*h, T2(i) + 0.5*k2_2); 
        k4_2 = h*F(Q2(i) + h, T2(i) + k3_2);
        T2(i+1) = T2(i) + (1/6)*(k1_2 + 2*k2_2 + 2*k3_2 + k4_2);
        %for Face 3
        k1_3 = h*F(Q3(i), T3(i)); 
        k2_3 = h*F(Q3(i) + 0.5*h, T3(i) + 0.5*k1_3); 
        k3_3 = h*F(Q3(i) + 0.5*h, T3(i) + 0.5*k2_3); 
        k4_3 = h*F(Q3(i) + h, T3(i) + k3_3);        
        T3(i+1) = T3(i) + (1/6)*(k1_3 + 2*k2_3 + 2*k3_3 + k4_3);
        %for Face 4
        k1_4 = h*F(Q4(i), T4(i)); 
        k2_4 = h*F(Q4(i) + 0.5*h, T4(i) + 0.5*k1_4); 
        k3_4 = h*F(Q4(i) + 0.5*h, T4(i) + 0.5*k2_4); 
        k4_4 = h*F(Q4(i) + h, T4(i) + k3_4);
        T4(i+1) = T4(i) + (1/6)*(k1_4 + 2*k2_4 + 2*k3_4 + k4_4);
        %for Face 5 
        k1_5 = h*F(Q5(i), T5(i)); 
        k2_5 = h*F(Q5(i) + 0.5*h, T5(i) + 0.5*k1_5); 
        k3_5 = h*F(Q5(i) + 0.5*h, T5(i) + 0.5*k2_5); 
        k4_5 = h*F(Q5(i) + h, T5(i) + k3_5);
        T5(i+1) = T5(i) + (1/6)*(k1_5 + 2*k2_5 + 2*k3_5 + k4_5);
        %for Face 6        
        k1_6 = h*F(Q6(i), T6(i)); 
        k2_6 = h*F(Q6(i) + 0.5*h , T6(i) + 0.5*k1_6); 
        k3_6 = h*F(Q6(i) + 0.5*h , T6(i) + 0.5*k2_6); 
        k4_6 = h*F(Q6(i) + h , T6(i) + k3_6);
        T6(i+1) = T6(i) + (1/6)*(k1_6 + 2*k2_6 + 2*k3_6 + k4_6);
 end
% t=1:7200;
% figure(1);
% legend('T6','T5','T4','T3','T2','T1');
% plot(T6,'r')
% hold on
% plot(T5,'g')
% plot(T4,'m')
% plot(T3,'y')
% plot(T2,'b')
% plot(T1,'k') 
% ylabel('Temperature in Spring');
% xlabel('Time in minutes');
