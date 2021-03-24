function [m] = magnetometer(qib,im)
%MAGNETOMETER Function to simulate the magnetometer on an IMU
%   Detailed explanation goes here

%Finding the rotational matrix from the quaternion 
qs = qib(1);
qx = qib(2);
qy = qib(3);
qz = qib(4);
%Using a loop to find R for every instance of qib
r11 = 1 - 2*qy^2 - 2*qz^2;
r12 = 2*qx*qy - 2*qz*qs;
r13 = 2*qx*qz + 2*qy*qs;
r21 = 2*qx*qy + 2*qz*qs;
r22 = 1 - 2*qx^2 - 2*qz^2;
r23 = 2*qy*qz - 2*qx*qs;
r31 = 2*qx*qz - 2*qy*qs  ;
r32 = 2*qy*qz + 2*qx*qs;
r33 = 1 - 2*qx^2 - 2*qy^2;

R = [r11 r12 r13;r21 r22 r23;r31 r32 r33]; % Rotation matrix

bm = R*im;  %Actual magnetic field in body frame

rms = 1000 * 1e-9;
rms_noise = randn(3,1)*rms;
% Init_cal_tol = 500;
% cross_talk = 0.6;

m = bm + rms_noise;
end

