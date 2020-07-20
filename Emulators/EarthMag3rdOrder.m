function [Bi] = EarthMag3rdOrder(r,t)
% EARTHMAG3RDORDER generates an EMF vector based upon a 3rd order model of
% the EMF adapted from Ryan Caverly's EarthMagField. Note that this is
% working code and is not necessarily complete.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('MagFieldCoeff_3rd_Order.mat')
g = g.*10^9;
h = h.*10^9;

a = 6371000; % Earth's Radius (m)
r1 = r(1);
r2 = r(2);
r3 = r(3);
rmag = sqrt(r1^2+r2^2+r3^2);
gNew = zeros(3,4); % Initialize new gnm values
hNew = zeros(3,4); % Initialize new hnm values

% Calculate right ascension of Greenwich Meridian
alphag0 = 0; % Initial RA of GM (rad)
alphag = alphag0 + 360.9856469*pi/180/86164*t + 0; % 23.9344696*3600 = 86164.0, use sidereal day
alphag = mod(alphag,2*pi); % Current RA of GM (rad)

% Calculate orbital angles of satellite
alpha = atan2(r2,r1); % Right ascension of satellite (rad)
delta = atan2(r3,sqrt(r1^2 + r2^2)); % Declination of satellite (rad)
phi = alpha - alphag; % Longitude of Satellite (rad)
theta = pi/2 - delta; % Coelevation of Satellite (rad)

sMat = SOrder3(); % Get S matrix
[pMat,dpMat] = POrder3(theta); % Get P and dP matrices

for i = 1:3 % n
    for j = 0:i % m
        gNew(i,j+1) = sMat(i+1,j+1)*g(i,j+1);
        hNew(i,j+1) = sMat(i+1,j+1)*h(i,j+1);
    end
end

Br = 0; % Initialize Br
for k = 1:3
    miniSum = 0; % Start Sum inside of sum
    for L = 0:k
        miniSum = miniSum + ((gNew(k,L+1)*cosd(L*phi) + hNew(k,L+1)*sind(L*phi))*pMat(k+1,L+1)); % Sum inside of sum
    end
    Br = Br + ((a/rmag)^(k+2)*(k+1)*miniSum); % Update Br
end

Btheta = 0; % Initialize Btheta
for k = 1:3
    miniSum = 0; % Start Sum inside of sum
    for L = 0:k
        miniSum = miniSum + ((gNew(k,L+1)*cosd(L*phi) + hNew(k,L+1)*sind(L*phi))*dpMat(k+1,L+1)); % Sum inside of sum
    end
    Btheta = Btheta - ((a/rmag)^(k+2)*miniSum); % Update Btheta
end

Bphi = 0; % Initialize Btheta
for k = 1:3
    miniSum = 0; % Start Sum inside of sum
    for L = 0:k
        miniSum = miniSum + L*((-gNew(k,L+1)*sind(L*phi) + hNew(k,L+1)*cosd(L*phi))*pMat(k+1,L+1)); % Sum inside of sum
    end
    Bphi = Bphi -((a/rmag)^(k+2)*miniSum); % Update Btheta
end
Bphi = (1/sind(theta))*Bphi; % Final update to Btheta
        
        
% Calculate ECI components of magnetic field vector
Bi1 = (Br*cos(delta) + Btheta*sin(delta))*cos(alpha) - Bphi*sin(alpha); %(T)
Bi2 = (Br*cos(delta) + Btheta*sin(delta))*sin(alpha) + Bphi*cos(alpha); %(T)
Bi3 = Br*sin(delta) - Btheta*cos(delta); %(T)

Bi = [Bi1 Bi2 Bi3]'; % Note, this is in the ECI frame!
        
end