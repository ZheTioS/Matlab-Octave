stpSize=0.1;
timeRun=[stpSize:stpSize:time.acc];
deltaOmega = alphaSat2*timeRun; % rad/sec
deltaSpeed = (rad2deg(deltaOmega)/360)*60; % rpm
newSpeed = sat2.pts.speed(end) + deltaSpeed;
a=length(sat2.pts.speed)+1;b=length(sat2.pts.speed)+length(newSpeed);
sat2.pts.speed(a:b)=newSpeed;
sat2.pts.time(a:b)=sat2.pts.time(end)+timeRun;
sat2.pts.angMom(a:b)=sat2.pts.angMom(end)+deltaOmega * sat2.MoI(1);