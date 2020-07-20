%% Quat Multiplication Test Script

% Define True 321 Euler Angles. Angles rotate from inertial frame to body frame.
Eul_b_I = [63 22 30]';          % [yaw pitch roll]', deg
Eul_b_Ir = Eul_b_I*pi/180;      % rad
C_b_I = eul2dcm(Eul_b_Ir);      % DCM from Inertial frame to body frame, True Angles (V_b = C_b_I*V_I) 
q_b_I = dcm2quat(C_b_I)';       % Quaternion from inertial frame to body frame

% Define Error 321 Euler angles.
Eul_Ihat_I = [1 2 3]';              % [yaw pitch roll]', deg
Eul_Ihat_Ir = Eul_Ihat_I*pi/180;    % rad
C_Ihat_I = eul2dcm(Eul_Ihat_Ir);    % DCM from Inertial frame to Estimate inertial frame
q_Ihat_I = dcm2quat(C_Ihat_I)';     % Quaternion from inertial frame to estimate inertial frame

% Calculate DCM from estimate inertial frame to body frame.
C_b_Ihat = C_b_I*C_Ihat_I';             % DCM from estimate inertial frame to body frame
q_b_Ihat = dcm2quat(C_b_Ihat)';         % Quaternion from estimate inertial frame to body frame
eul_b_Ihat=dcm2eul(quat2dcm(q_b_Ihat'));% Euler angles corresponding to rotation from estimate inertial frame to body frame
eul_b_Ihat=eul_b_Ihat*180/pi            % deg

% Test Quaternion Multiplication          
q_b_Ihat1 = quatmult(quatInverse(q_Ihat_I'),q_b_I); % Quaternion using multiplication orientation #1
eul_b_Ihat1=dcm2eul(quat2dcm(q_b_Ihat1'));          % Euler angles from quaternion
eul_b_Ihat1=eul_b_Ihat1*180/pi                      % deg
q_b_Ihat2 = quatmult(q_b_I,quatInverse(q_Ihat_I')); % Quaternion using multiplication orientation #2
eul_b_Ihat2=dcm2eul(quat2dcm(q_b_Ihat2'));          % Euler angles from quaternion 
eul_b_Ihat2=eul_b_Ihat2*180/pi                      % deg

%%
% We see that the quaternion multiplication method of q_b_Ihat2 gives the
% same value as our DCM calculation given by C_b_Ihat. This implies that
% the quaternion multiplication is in the same order as DCM convention.
%% Quat Error Test
q_Ihat_IT1 = quatmult(q_b_I,quatInverse(q_b_Ihat2'))
q_Ihat_IT2 = quatmult(quatInverse(q_b_Ihat2'),q_b_I)
q_Ihat_I
