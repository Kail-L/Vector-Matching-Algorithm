function [sMat] = SOrder3()

% [sMat] = SOrder3()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates Snm values for up to n = 3
%
% Inputs:
% None
%
% Outputs:
% sMat = 4x4 S matrix
%      NOTE: sMat(i,i) = S(i-1,i-1), for example: sMat(1,1) = S00, sMat(2,1) = S10, etc.
% 
% Joel Luedke
% Updated 04/15/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize S matrix
sMat = zeros(4,4);
sMat(0+1,0+1) = 1; % S00, doesn't depend on other values of S

for i=1:3 % n
    for j = 0:i % m
        if j == 0
           sMat(i+1,j+1) = sMat(i,j+1)*((2*i-1)/i); % For SX0
        else
           if j==1
               kDelta1 = 1; % Kronecker delta, 1 if m==1
           else
               kDelta1 = 0; % 0 if m/=1
           end
           sMat(i+1,j+1) = sMat(i+1,j)*sqrt((i-j+1)*(kDelta1+1)/(i+j)); % For Snm
        end
        
    end
end
end