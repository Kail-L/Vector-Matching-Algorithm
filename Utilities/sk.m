 function as=sk(a)
%SK  Computes skew symmetric matrix for 3x1 vector.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AS = sk(a) This function puts a 3x1 vector into skey symmetric format.
% 
% SOURCES:
% n/a (see skey symmetric matrices in google search)
%
% INPUT PARAMETERS:
% a = 3x1 vector
%
% OUTPUT PARAMETERS:
% as = 3x3 skey symmetric matrix
%
% VARIABLES:
% n/a
%
% Kail Laughlin
% Updated 11/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

as=[0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];

 end
