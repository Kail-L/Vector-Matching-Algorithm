function [pMat,dpMat] = POrder3(theta)

% [pMat,dpMat] = POrder3(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates P and dP matrices for a given value of theta. Used in
% EarthMag3rdOrder.m
%
% Inputs:
% theta = Coelevation angle of spacecraft (radians)
%
% Outputs:
% pMat = 4x4 P matrix
% dPmat = 4x4 dP matrix
%       NOTE: pMat(i,i) = P(i-1,i-1), for example: pMat(1,1) = P00, pMat(2,1) = P10, etc.
%       Same for dP matrix
%
% Joel Luedke
% Updated 04/15/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% Initialize P matrix
pMat = zeros(4,4);
pMat(0+1,0+1) = 1; % P00, doesn't depend on other P values
pMat(1+1,0+1) = cos(theta); % P10, doesn't depend on other P values
pMat(1+1,1+1) = sin(theta); % P11, doesn't depend on other P values

% Initialize dP matrix
dpMat = zeros(4,4);
dpMat(0+1,0+1) = 0; % dP00, doesn't depend on other dP values
dpMat(1+1,0+1) = -sin(theta)*pMat(0+1,0+1); % dP10
dpMat(1+1,1+1) = cos(theta)*pMat(0+1,0+1); % dP11


% Initialize K matrix
kMat = zeros(4,4);

for i = 2:3 % n
    for j = 0:i-1 % m
    % Fill in K, P, and dP matrices for n/=m
    kMat(i+1,j+1) = ((i-1)^2-(j^2))/((2*i-1)*(2*i-3));
    pMat(i+1,j+1) = cos(theta)*pMat(i,j+1) - kMat(i+1,j+1)*pMat(i-1,j+1);
    dpMat(i+1,j+1) = cos(theta)*dpMat(i,j+1) - sin(theta)*pMat(i,j+1) -...
        kMat(i+1,j+1)*dpMat(i-1,j+1);
    
    % Fill in dP
    end
    % Fill in P and dP matrices for n==m
    pMat(i+1,i+1) = sin(theta)*pMat(i,i);
    dpMat(i+1,i+1) = sin(theta)*dpMat(i,i) + cos(theta)*pMat(i,i);
end
end