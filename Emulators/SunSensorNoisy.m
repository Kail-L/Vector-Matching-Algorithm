function [vector_noisy,vector_clean] = SunSensorNoisy(vector,Cba,const)
%SUNSENSORNOISY   Calculate noisy sun sensor measurement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [vector_noisy,vector_clean] = SunSensorNoisy(vector,Cba,t) determines a noisy 
% measurement from a sun sensor while also providing a clean measurement as
% well.
%
% SOURCES:
% http://folk.ntnu.no/tomgra/Diplomer/Sunde.pdf, pg. 29
%
% INPUT PARAMETERS:
% vector = 3x1 vector in a frame
% Cba = 3x3 DCM of frame b relative to frame a
% t = time
%
% OUTPUT PARAMETERS:
% vector_noisy = 3x1 noisy unit vector measurement in b frame
%
% Kail Laughlin
% Updated 1/28/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% error in each angle
err=const.S_measNoise;
err_angle = err*pi/180; % rad

% noise in pitch, roll, and yaw
noise_pitch = err_angle*(2*rand-1);
noise_roll = err_angle*(2*rand-1);
noise_yaw = err_angle*(2*rand-1);

% corrupting rotation matrix
EulCorrup=[noise_yaw;
           noise_pitch;
           noise_roll];

% error rotation matrix
Cbae=eul2dcm(EulCorrup);

if norm(vector)==0
    vector_noisy = EulCorrup;
    vector_clean = Cba*vector';
else
    vector_noisy = Cbae*Cba*vector';
    vector_clean = Cba*vector';
end
end