function Gq=QuatJoc(q)
%QUATJOC Function to calculate quaternion jacobian w.r.t. Euler Angles.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [GQ] = QUATJOC(Q) This function uses the quaternion in the form of
% [scalar; vector] to determine the jacobian matrix with respect to the
% euler angles given by the quaternion. This matrix is used to convert the
% covariance on the euler angles to a covariance on the euler angles for a
% given step through the following equation:
%
% P_eul = Gq*P_quat*Gq'
%
% SOURCES:
% Main equation taken from the Thesis of John B. Schleppe (see address)
% www.ucalgary.ca/engo_webdocs/GL/96.20096.JSchleppe.pdf Pg. 69
%
% Our Quaternion derivatives taken from Wikipedia:
% https://en.wikipedia.org/wiki/Conversion_between_quaternions_and_Euler_angles
%
% Additionally, partial derivatives were calculated using the following:
% 
% syms q0 q1 q2 q3
% assume(q0,'real')
% assume(q1,'real')
% assume(q2,'real')
% assume(q3,'real')
% phi=atan2(2*(q0*q1+q2*q3),1-2*(q1^2+q2^2));
% the=asin(2*(q0*q2-q3*q1));
% psi=atan2(2*(q0*q3+q1*q2),1-2*(q2^2+q3^2));
% dpsi_dq0 = diff(psi,q0);
% dpsi_dq1 = diff(psi,q1);
% dpsi_dq2 = diff(psi,q2);
% dpsi_dq3 = diff(psi,q3);
% dtheta_dq0 = diff(the,q0);
% dtheta_dq1 = diff(the,q1);
% dtheta_dq2 = diff(the,q2);
% dtheta_dq3 = diff(the,q3);
% dphi_dq0 = diff(phi,q0);
% dphi_dq1 = diff(phi,q1);
% dphi_dq2 = diff(phi,q2);
% dphi_dq3 = diff(phi,q3);
% 
% G = [dpsi_dq0 dpsi_dq1 dpsi_dq2 dpsi_dq3;
%      dtheta_dq0 dtheta_dq1 dtheta_dq2 dtheta_dq3;
%      dphi_dq0 dphi_dq1 dphi_dq2 dphi_dq3];
%
% INPUT PARAMETERS:
% q = 1x4 quaternion with the scalar component first in form of [scalar;
% vector].
%
% OUTPUT PARAMETERS:
% Gq = 3x4 jacobian matrix of euler angles with respect to quaternion.
%
% VARIABLES:
% q0 = first value of quaternion (scalar).
% q1 = second value of quaternion (vector).
% q2 = third value of quaternion (vector).
% q3 = fourth value of quaternion (vector).
% dpsi_dq0 = partial derivative of psi w.r.t q1.
% dpsi_dq1 = partial derivative of psi w.r.t q2.
% dpsi_dq2 = partial derivative of psi w.r.t q3.
% dpsi_dq3 = partial derivative of psi w.r.t q4.
% dtheta_dq0 = partial derivative of theta w.r.t q1.
% dtheta_dq1 = partial derivative of theta w.r.t q2.
% dtheta_dq2 = partial derivative of theta w.r.t q3.
% dtheta_dq3 = partial derivative of theta w.r.t q4.
% dphi_dq0 = partial derivative of phi w.r.t q1.
% dphi_dq1 = partial derivative of phi w.r.t q2.
% dphi_dq2 = partial derivative of phi w.r.t q3.
% dphi_dq3 = partial derivative of phi w.r.t q4.
%
% Created by: Kail Laughlin
% Updated by: Kail Laughlin
% Last Modified Date: 3/6/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q0=q(1);
q1=q(2);
q2=q(3);
q3=q(4);

dpsi_dq0 = -(2*q3*(2*q2^2 + 2*q3^2 - 1))/((2*q2^2 + 2*q3^2 - 1)^2 + (2*q0*q3 + 2*q1*q2)^2);
dpsi_dq1 = -(2*q2*(2*q2^2 + 2*q3^2 - 1))/((2*q2^2 + 2*q3^2 - 1)^2 + (2*q0*q3 + 2*q1*q2)^2);
dpsi_dq2 = -(((2*q1)/(2*q2^2 + 2*q3^2 - 1) - (4*q2*(2*q0*q3 + 2*q1*q2))/(2*q2^2 + 2*q3^2 - 1)^2)*(2*q2^2 + 2*q3^2 - 1)^2)/((2*q2^2 + 2*q3^2 - 1)^2 + (2*q0*q3 + 2*q1*q2)^2);
dpsi_dq3 = -(((2*q0)/(2*q2^2 + 2*q3^2 - 1) - (4*q3*(2*q0*q3 + 2*q1*q2))/(2*q2^2 + 2*q3^2 - 1)^2)*(2*q2^2 + 2*q3^2 - 1)^2)/((2*q2^2 + 2*q3^2 - 1)^2 + (2*q0*q3 + 2*q1*q2)^2);

dtheta_dq0 = (2*q2)/(1 - (2*q0*q2 - 2*q1*q3)^2)^(1/2);
dtheta_dq1 = -(2*q3)/(1 - (2*q0*q2 - 2*q1*q3)^2)^(1/2);
dtheta_dq2 = (2*q0)/(1 - (2*q0*q2 - 2*q1*q3)^2)^(1/2);
dtheta_dq3 = -(2*q1)/(1 - (2*q0*q2 - 2*q1*q3)^2)^(1/2);

dphi_dq0 = -(2*q1*(2*q1^2 + 2*q2^2 - 1))/((2*q1^2 + 2*q2^2 - 1)^2 + (2*q0*q1 + 2*q2*q3)^2);
dphi_dq1 = -(((2*q0)/(2*q1^2 + 2*q2^2 - 1) - (4*q1*(2*q0*q1 + 2*q2*q3))/(2*q1^2 + 2*q2^2 - 1)^2)*(2*q1^2 + 2*q2^2 - 1)^2)/((2*q1^2 + 2*q2^2 - 1)^2 + (2*q0*q1 + 2*q2*q3)^2);
dphi_dq2 = -(((2*q3)/(2*q1^2 + 2*q2^2 - 1) - (4*q2*(2*q0*q1 + 2*q2*q3))/(2*q1^2 + 2*q2^2 - 1)^2)*(2*q1^2 + 2*q2^2 - 1)^2)/((2*q1^2 + 2*q2^2 - 1)^2 + (2*q0*q1 + 2*q2*q3)^2);
dphi_dq3 = -(2*q2*(2*q1^2 + 2*q2^2 - 1))/((2*q1^2 + 2*q2^2 - 1)^2 + (2*q0*q1 + 2*q2*q3)^2);

% Gq = [dpsi_dq0 dpsi_dq1 dpsi_dq2 dpsi_dq3;
%       dtheta_dq0 dtheta_dq1 dtheta_dq2 dtheta_dq3;
%       dphi_dq0 dphi_dq1 dphi_dq2 dphi_dq3];
  
Gq = [dphi_dq0 dphi_dq1 dphi_dq2 dphi_dq3;      
      dtheta_dq0 dtheta_dq1 dtheta_dq2 dtheta_dq3;
      dpsi_dq0 dpsi_dq1 dpsi_dq2 dpsi_dq3];
end
