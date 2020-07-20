function Eul = quat2eul(q);
%QUAT2EUL Function to calculate 3-2-1 Euler Angles from Quaternion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eul = QUAT2EUL(q) This function uses the quaternion in the form of
% [scalar; vector] to determine the 3-2-1 Euler angles for the given
% quaternion rotation. We utilize the equations given below:
%
% [psi theta phi]^T = [arctan(2(q0q3+q1q2)/(1-2(q2^2+q3^2)));
%                      arcsin(2(q0q2-q3q1));
%                      arctan(2(q0q1+q2q3)/(1-2(q1^2+q2^2)))]
%
% Note: When implementing in computer format, we use atan2:
%
% [psi theta phi]^T = [atan2(2(q0q3+q1q2),1-2(q2^2+q3^2));
%                      asin(2(q0q2-q3q1));
%                      atan2(2(q0q1+q2q3),1-2(q1^2+q2^2))]
% SOURCES:
% Wikepedia (Its correct)
% https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
%
% INPUT PARAMETERS:
% q = 4x1 quaternion with the scalar component first in form of [scalar;
% vector].
%
% OUTPUT PARAMETERS:
% Eul = 3x1 array of euler angles (rad).
%
% VARIABLES:
% q0 = first value of quaternion (scalar).
% q1 = second value of quaternion (vector).
% q2 = third value of quaternion (vector).
% q3 = fourth value of quaternion (vector).
%
% Created by: Kail Laughlin
% Updated by: Kail Laughlin
% Last Modified Date: 4/20/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q0 = q(1);
q1 = q(2);
q2 = q(3);
q3 = q(4);

Eul = [atan2(2*(q0*q3+q1*q2),1-2*(q2^2+q3^2));
       asin(2*(q0*q2-q3*q1));
       atan2(2*(q0*q1+q2*q3),1-2*(q1^2+q2^2))];
   
end
   