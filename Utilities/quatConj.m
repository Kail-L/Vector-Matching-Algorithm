function conj = quatConj(q)
%QUATCONJ function to determine conjugate of quaternion.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONJ = quatConj(q) returns the conjugate of an input quaternion.
%
% INPUT PARAMETERS:
% q = 4x1 quaternion
%
% OUTPUT PARAMETERS:
% conj = 4x1 conjugate quaternion
%
% Kail Laughlin
% Updated 11/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
conj = [q(1); -q(2:4)];

end