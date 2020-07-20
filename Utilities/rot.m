function R = rot(theta,primary)
%ROT  Computes the primary rotation matrix given an angle of rotation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% R = rot(theta,primary) This function creates a rotation matrix for a
% single primary rotation (1 rotation, 2 rotation, or 3 rotation) given an
% angle and the axis of rotation.
% 
% SOURCES:
% n/a (see primary rotation matrices in google search)
%
% INPUT PARAMETERS:
% theta = scalar angle of rotation (rad)
% primary = integer specifying axis of rotation
% NOTE: primary must be scalar 1, 2, or 3
%
% OUTPUT PARAMETERS:
% R = 3x3 rotation matrix corresponding to angle and axis given as inputs
%
% VARIABLES:
% n/a
%
% Kail Laughlin
% Updated 11/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch primary
    case 1
        R = [ 1 0 0 ; 0 cos(theta) sin(theta); 0 -sin(theta) cos(theta)];
    case 2 
        R = [ cos(theta) 0 -sin(theta) ; 0 1 0; sin(theta) 0 cos(theta)];
    case 3
        R = [cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0 ; 0 0 1];
end