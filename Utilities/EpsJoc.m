function Eq=EpsJoc(q)
%EPSJOC Function to calculate quaternion jacobian w.r.t. Epsilon.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [EQ] = QUATJOC(Q) This function uses the quaternion in the form of
% [scalar; vector] to determine the jacobian matrix with respect to
% epsilon given by the quaternion. This matrix is used to convert the
% covariance on the quaternion to a covariance on epsilon for a
% given step through the following equation:
%
% P_eps = Eq*P_quat*Eq'
%
% SOURCES:
% Main equation adapted from the Thesis of John B. Schleppe (see address)
% www.ucalgary.ca/engo_webdocs/GL/96.20096.JSchleppe.pdf Pg. 69
%
% The equation to calculate epsilon is given by
%
%       Epsilon = acos(Cbbhat(3,3)/(norm(Cbbhat(:,3)))*180/pi;
%
% where Cbbhat is the error dcm created based on the error quaternion.
% Specifically, 
%
%       Cbbhat(3,3) = dq0^2-dq1^2-dq2^2-dq3^2;
%       Cbbhat(:,3) = [2(dq1dq3-dq0dq2);
%                      2(dq2dq3+dq0dq1);
%                      dq0^2-dq1^2-dq2^2-dq3^2]
%
% with dq0, dq1, dq2, & dq3 as the delta quaternions calculated from the
% EKF. Since we have a covariance on these states, we just need to use the
% P_eps equation given above to convert to a covariance on epsilon.
%
% Additionally, partial derivatives were calculated using the following:
% 
% syms q0 q1 q2 q3
% assume(q0,'real')
% assume(q1,'real')
% assume(q2,'real')
% assume(q3,'real')
% 
% Cbbhat33 = q0^2-q1^2-q2^2+q3^2;
% Cbbhat3v = [2*(q1*q3-q0*q2);
%             2*(q2*q3+q0*q1);
%             q0^2-q1^2-q2^2+q3^2];
%         
% epsilon = acos(Cbbhat33/norm(Cbbhat3v));
% deps_dq0 = diff(epsilon,q0);
% deps_dq1 = diff(epsilon,q1);
% deps_dq2 = diff(epsilon,q2);
% deps_dq3 = diff(epsilon,q3);
% 
% INPUT PARAMETERS:
% q = 1x4 quaternion with the scalar component first in form of [scalar;
% vector].
%
% OUTPUT PARAMETERS:
% Eq = 1x4 jacobian matrix of euler angles with respect to quaternion.
%
% VARIABLES:
% q0 = first value of quaternion (scalar).
% q1 = second value of quaternion (vector).
% q2 = third value of quaternion (vector).
% q3 = fourth value of quaternion (vector).
% deps_dq0 = partial derivative of epsilon w.r.t q1.
% deps_dq1 = partial derivative of epsilon w.r.t q2.
% deps_dq2 = partial derivative of epsilon w.r.t q3.
% deps_dq3 = partial derivative of epsilon w.r.t q4.
%
% Created by: Kail Laughlin
% Updated by: Kail Laughlin
% Last Modified Date: 4/20/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q0=q(1);
q1=q(2);
q2=q(3);
q3=q(4);



deps_dq0 =-((2*q0)/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(1/2) - ((4*q0*abs(q0^2 - q1^2 - q2^2 + q3^2)*sign(q0^2 - q1^2 - q2^2 + q3^2) + 4*q1*abs(2*q0*q1 + 2*q2*q3)*sign(2*q0*q1 + 2*q2*q3) + 4*q2*abs(2*q0*q2 - 2*q1*q3)*sign(2*q0*q2 - 2*q1*q3))*(q0^2 - q1^2 - q2^2 + q3^2))/(2*(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(3/2)))/(1 - (q0^2 - q1^2 - q2^2 + q3^2)^2/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2))^(1/2);
deps_dq1 =((2*q1)/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(1/2) - ((4*q1*abs(q0^2 - q1^2 - q2^2 + q3^2)*sign(q0^2 - q1^2 - q2^2 + q3^2) - 4*q0*abs(2*q0*q1 + 2*q2*q3)*sign(2*q0*q1 + 2*q2*q3) + 4*q3*abs(2*q0*q2 - 2*q1*q3)*sign(2*q0*q2 - 2*q1*q3))*(q0^2 - q1^2 - q2^2 + q3^2))/(2*(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(3/2)))/(1 - (q0^2 - q1^2 - q2^2 + q3^2)^2/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2))^(1/2);
deps_dq2 = ((2*q2)/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(1/2) + ((4*q0*abs(2*q0*q2 - 2*q1*q3)*sign(2*q0*q2 - 2*q1*q3) - 4*q2*abs(q0^2 - q1^2 - q2^2 + q3^2)*sign(q0^2 - q1^2 - q2^2 + q3^2) + 4*q3*abs(2*q0*q1 + 2*q2*q3)*sign(2*q0*q1 + 2*q2*q3))*(q0^2 - q1^2 - q2^2 + q3^2))/(2*(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(3/2)))/(1 - (q0^2 - q1^2 - q2^2 + q3^2)^2/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2))^(1/2);
deps_dq3 = -((2*q3)/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(1/2) - ((4*q3*abs(q0^2 - q1^2 - q2^2 + q3^2)*sign(q0^2 - q1^2 - q2^2 + q3^2) - 4*q1*abs(2*q0*q2 - 2*q1*q3)*sign(2*q0*q2 - 2*q1*q3) + 4*q2*abs(2*q0*q1 + 2*q2*q3)*sign(2*q0*q1 + 2*q2*q3))*(q0^2 - q1^2 - q2^2 + q3^2))/(2*(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(3/2)))/(1 - (q0^2 - q1^2 - q2^2 + q3^2)^2/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2))^(1/2);

% deps_dq0 =-(180*((2*q0)/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(1/2) - ((4*q0*abs(q0^2 - q1^2 - q2^2 + q3^2)*sign(q0^2 - q1^2 - q2^2 + q3^2) + 4*q1*abs(2*q0*q1 + 2*q2*q3)*sign(2*q0*q1 + 2*q2*q3) + 4*q2*abs(2*q0*q2 - 2*q1*q3)*sign(2*q0*q2 - 2*q1*q3))*(q0^2 - q1^2 - q2^2 + q3^2))/(2*(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(3/2))))/(pi*(1 - (q0^2 - q1^2 - q2^2 + q3^2)^2/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2))^(1/2));
% deps_dq1 =(180*((2*q1)/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(1/2) - ((4*q1*abs(q0^2 - q1^2 - q2^2 + q3^2)*sign(q0^2 - q1^2 - q2^2 + q3^2) - 4*q0*abs(2*q0*q1 + 2*q2*q3)*sign(2*q0*q1 + 2*q2*q3) + 4*q3*abs(2*q0*q2 - 2*q1*q3)*sign(2*q0*q2 - 2*q1*q3))*(q0^2 - q1^2 - q2^2 + q3^2))/(2*(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(3/2))))/(pi*(1 - (q0^2 - q1^2 - q2^2 + q3^2)^2/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2))^(1/2));
% deps_dq2 =(180*((2*q2)/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(1/2) + ((4*q0*abs(2*q0*q2 - 2*q1*q3)*sign(2*q0*q2 - 2*q1*q3) - 4*q2*abs(q0^2 - q1^2 - q2^2 + q3^2)*sign(q0^2 - q1^2 - q2^2 + q3^2) + 4*q3*abs(2*q0*q1 + 2*q2*q3)*sign(2*q0*q1 + 2*q2*q3))*(q0^2 - q1^2 - q2^2 + q3^2))/(2*(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(3/2))))/(pi*(1 - (q0^2 - q1^2 - q2^2 + q3^2)^2/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2))^(1/2));
% deps_dq3 =-(180*((2*q3)/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(1/2) - ((4*q3*abs(q0^2 - q1^2 - q2^2 + q3^2)*sign(q0^2 - q1^2 - q2^2 + q3^2) - 4*q1*abs(2*q0*q2 - 2*q1*q3)*sign(2*q0*q2 - 2*q1*q3) + 4*q2*abs(2*q0*q1 + 2*q2*q3)*sign(2*q0*q1 + 2*q2*q3))*(q0^2 - q1^2 - q2^2 + q3^2))/(2*(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2)^(3/2))))/(pi*(1 - (q0^2 - q1^2 - q2^2 + q3^2)^2/(abs(2*q0*q1 + 2*q2*q3)^2 + abs(2*q0*q2 - 2*q1*q3)^2 + abs(q0^2 - q1^2 - q2^2 + q3^2)^2))^(1/2));

% deps_dq0 =0;
% deps_dq1 =(2*q1)/(1 - (q1^2 + q2^2 - q3^2 - 1)^2)^(1/2);
% deps_dq2 =(2*q2)/(1 - (q1^2 + q2^2 - q3^2 - 1)^2)^(1/2);
% deps_dq3 =-(2*q3)/(1 - (q1^2 + q2^2 - q3^2 - 1)^2)^(1/2);

Eq = [deps_dq0 deps_dq1 deps_dq2 deps_dq3];

end
