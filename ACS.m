% w=zeros(1,11);
% t=0:10;
% I=7*10^(-7);
% w(1)=12000;
% alpha = 1882.7;
% T = 23*10^(-6);
% A = T/I;
% q(1)=T;
% q(2)=T;
% q(3)=T;
% q(4)=T;
% q(5)=T;
% q(6)=T;
% q(7)=-T;
% q(8)=-T;
% q(9)=-T;
% q(10)=-T;
% q(11)=-T;
% 
% w(2)=(w(1)+(alpha*1));
% w(3)=(w(1)+(alpha*2));
% w(4)=(w(1)+(alpha*3));
% w(5)=(w(1)+(alpha*4));
% w(6)=(w(1)+(alpha*5));
% w(7)=(w(6)-(alpha*1));
% w(8)=(w(6)-(alpha*2));
% w(9)=(w(6)-(alpha*3));
% w(10)=(w(6)-(alpha*4));
% w(11)=(w(6)-(alpha*5));
% 
% for i=1:11
% H(i)=I*w(i);
% end
% % plot(w)
% % ylabel('Angular velocity')
% fig1=figure;
% hold on
% % plot(t,q)
%  plot(t,H)
% ylabel('Angular Momentum')
% % plot([0 5 5 10], T*[-1 -1 1 1])
% plot(t,-H)
% xlabel('time')
% plot(t,H+(-H))
% legend('Reaction Wheel','Spacecraft','System','location','best')
% hold off
% fig2=figure
% Sat.mass    =1  ;         %   kg
% Sat.Length  =0.1;         %   m
% Sat.Height  =0.1;         %   m
% Sat.width   =0.1;         %   m
% Sat.torque=23*power(10,-6);
% Sat.I=7*power(10,-7);
% 
% Sat.area=Sat.Height*Sat.Length;    %   m^2
% 
% Sat.MOI=(1/6)*Sat.mass*Sat.area;
% 
% Sat.alpha=(Sat.torque/Sat.MOI);
% 
% Sat.N=Sat.MOI*Sat.alpha;
% Sat.accl=Sat.N/Sat.MOI;
% t=[1:10];
% Sat.theta=0.5*(Sat.accl)*power(t,2);
% Sat.angle_In_Rads=((180./pi).*((t.^2).*Sat.torque)/(Sat.MOI.*2)).*ones(1, length(t)) ;  % Angle in rads
% plot(t,Sat.angle_In_Rads)
% title('Maneuver Angle');
% ylabel('Angle');
% xlabel('Time');
% % grid on
% % fig3=figure;
% % hold on
% % plot(t, H)
% % ylabel('Angular Momentum of Reaction Wheel')
% % plot(t, -H)
% % plot(sat2.pts.time, sat2.pts.angMom+sat2.rw.pts.angMom)
% % % ylabel('Spacecraft Angular Momentum, kg-m2/sec')
% % legend('Reaction Wheel','Spacecraft','location','best')
% % xlabel('Time')
% % hold off
% % % x = w(1)+(alpha*5);

alphaRw = (4*10^(-4))/(1.2*10^(-7)) ; %rad/sec^2

% Calculate the time it takes to get to nomimal speed
newSpeed = 8000 - sat3.rw.pts.speed(end)  ;
newOmega=deg2rad(newSpeed/60*360) ;
timeNeeded=newOmega/alphaRw ;

stpSize    = 0.01 ;                % sec
time.acc   = round(timeNeeded,2) ; % sec
time.dec   = time.acc ;            % sec
time.coast = 100 - 2*time.acc ;    % sec

alphaSat3 = sat3.rw.info.torqueNom/sat3.MoI(1);

% Accel
timeRun=[stpSize:stpSize:time.acc];
deltaOmega = alphaSat3*timeRun; % rad/sec
deltaSpeed = (rad2deg(deltaOmega)/360)*60; % rpm
newSpeed = sat3.pts.speed(end) + deltaSpeed;
a=length(sat3.pts.speed)+1;b=length(sat3.pts.speed)+length(newSpeed);
sat3.pts.speed(a:b)=newSpeed;
sat3.pts.time(a:b)=sat3.pts.time(end)+timeRun;
sat3.pts.angMom(a:b)=sat3.pts.angMom(end)+...
  deltaOmega * sat3.MoI(1);

% Hold Max Speed
timeRun=[stpSize:stpSize:time.coast];
newSpeed = sat3.pts.speed(end) * ones(length(timeRun),1) ;
a=length(sat3.pts.speed)+1;b=length(sat3.pts.speed)+length(newSpeed);
sat3.pts.speed(a:b)=newSpeed ;
sat3.pts.time(a:b)=sat3.pts.time(end)+timeRun;
sat3.pts.angMom(a:b)=sat3.pts.angMom(end)* ones(length(timeRun),1);

% Decel
timeRun=[stpSize:stpSize:time.dec];
deltaOmega = -alphaSat3*timeRun; % rad/sec
deltaSpeed = (rad2deg(deltaOmega)/360)*60; % rpm
newSpeed = sat3.pts.speed(end) + deltaSpeed; % rpm
a=length(sat3.pts.speed)+1;b=length(sat3.pts.speed)+length(newSpeed);
sat3.pts.speed(a:b)=newSpeed;
sat3.pts.time(a:b)=sat3.pts.time(end)+timeRun;
sat3.pts.angMom(a:b)=sat3.pts.angMom(end)+...
  deltaOmega * sat3.MoI(1);
