 function eul = dcm2eul(C)
%DCM2EUL  Converts DCM into 321 Euler angles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EUL = dcm2eul(C) This function converts a given DCM into 321 euler angles
% which rotate two given frames. Be careful with the notation of DCMs that
% you put in to ensure that the angles which come out are what are
% expected.
% 
% SOURCES:
% n/a (see any DCM to euler angle code)
%
% INPUT PARAMETERS:
% C = 3x3 direction cosine matrix (unitless)
%
% OUTPUT PARAMETERS:
% eul = 3x1 vector containing the yaw, pitch roll angles (radians)
%     = [psi theta phi]^T
%
% VARIABLES:
% n/a
%
% Kail Laughlin, taken from Demoz Gebre
% Updated 11/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eul(1,1) = atan2(C(1,2),C(1,1));  %     yaw
eul(2,1) = asin(-C(1,3));         %     pitch
eul(3,1) = atan2(C(2,3),C(3,3));  %     roll
 end
