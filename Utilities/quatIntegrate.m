function [q_gyro] = quatIntegrate(gyro_vec, q_prev, dt)
%QUATINTEGRATE function to integrate quaternion through time using gyro.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [q_gyro] = quatIntegrate(gyro_vec, q_prev, dt) returns the updated
% quaternion at the next time step (t_next = t_current + dt).
%
% SOURCES:
% Quaternion kinematics for the error-state Kalman filter, Joan Sola 2017
%
% INPUT PARAMETERS:
% gyro_vec = 3x1 gyro rates in body frame (rad/s)
% q_prev = 4x1 quaternion state
% dt = time step to integrate with (s)
%
% OUTPUT PARAMETERS:
% q_gyro = 4x1 updated quaternion using gyroscope readings.
%
% VARIABLES:
% gyro_skew_sym = 4x4 skew symmetric matrix using gyro values
% dq = 4x1 rate of change matrix of quaternion
%
% Kail Laughlin
% Updated 11/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Using Euler Step %
% 
% % create skew Symmetric matrix
% gyro_skew_sym = [0 -gyro_vec(1) -gyro_vec(2) -gyro_vec(3);
%                  gyro_vec(1) 0 gyro_vec(3) -gyro_vec(2);
%                  gyro_vec(2) -gyro_vec(3) 0 gyro_vec(1);
%                  gyro_vec(3) gyro_vec(2) -gyro_vec(1) 0];
%              
% %find ROC of q                
% dqdt = 1/2 * gyro_skew_sym * q_prev;
% 
% %find q_new
% q_gyro = q_prev+dqdt.*dt;

% Using Quaternion Multiplication %
dq = [cos(norm(gyro_vec)*dt/2);...
     (gyro_vec/norm(gyro_vec)*sin(norm(gyro_vec)*dt/2))'];

q_gyro = quatmult(dq,q_prev);

end