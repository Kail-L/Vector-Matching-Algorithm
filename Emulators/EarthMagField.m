function bi = EarthMagField(r,t)
%EARTHMAGFIELD  Calculates magnetic field vector at given point in orbit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BI = EarthMagField(r,t) This function computes the Earth's magnetic field
% physical vector resolved in the ECI frame given position vector of
% satellite in ECI frame and time.
% 
% SOURCES:
% James Richard Forbes, Ryan Caverly
% Also, Appendix H of Wertz (1978), Spacecraft Attitdue Determination and
% Control, Kluwer, Academic Publishers, Dordrecht, The Netherlands.
% https://link.springer.com/content/pdf/bbm%3A978-94-009-9907-7%2F1.pdf
% Page 779 (Pg 55 of the PDF)
%
% INPUT PARAMETERS:
% r = 3x1 position vector of satellite in ECI frame (m)
% t = current time in orbit (s)
%
% OUTPUT PARAMETERS:
% bi = 3x1 magnetic field vector at given point in orbit, ECI frame (T)
%
% VARIABLES:
% Re = radius of earth (m)
% g10 = magnetic field constants (T)
% g11 = magnetic field constants (T)
% h11 = magnetic field constants (T)
% r1 = "x" position of spacecraft in ECI frame (m)
% r2 = "y" position of spacecraft in ECI frame (m)
% r3 = "z" position of spacecraft in ECI frame (m)
% r_mag = magnitude of position vector (m)
% alphag0 = Initial Right Ascension of Greenwich Meridian (rad)
% alphag = Right Ascension of Greenwich Meridian (rad)
% alpha = Right Ascension of Satellite (rad)
% delta = Declination of Satellite (rad)
% phi = East Longitude from Greenwich of Satellite (rad)
% theta = Coelevation of Satellite (rad)
% br = Radial Component of Magnetic Field Vector (Outward Positive) (T)
% btheta = Coelevation Component of Magnetic Field Vector (South Positive) (T)
% bphi = Azimuthal Component of Magnetic Field Vector (East Positive) (t)
% bi1 = mag field component in r1 direction (T)
% bi2 = mag field component in r2 direction (T)
% bi3 = mag field component in r3 direction (T)
%
% Kail Laughlin, Joel Luedke, taken from James Richard Forbes
% Updated 02/13/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Re = 6371.2*1000; % m, Radius of the earth

% IGRF Magnetic field constants (2020) %https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
g10 = -29405e-9; % T 
g11 = -1451e-9; % T
h11 = 4653e-9; % T

% Position of satellite in ECI frame. 
r1 = r(1); % m
r2 = r(2); % m
r3 = r(3); % m
r_mag = sqrt(r.'*r); % m

% Calculate right ascension of Greenwich Meridian
alphag0 = 0; % Initial RA of GM (rad)
alphag = alphag0 + 360.9856469*pi/180/86164*t + 0; % 23.9344696*3600 = 86164.0, use sidereal day
alphag = mod(alphag,2*pi); % Current RA of GM (rad)

% Calculate orbital angles of satellite
alpha = atan2(r2,r1); % Right ascension of satellite (rad)
delta = atan2(r3,sqrt(r1^2 + r2^2)); % Declination of satellite (rad)
phi = alpha - alphag; % Longitude of Satellite (rad)
theta = pi/2 - delta; % Coelevation of Satellite (rad)

% Calculate radial, coelevation, and azimuthal components of magnetic field vector
br = 2*(Re/r_mag)^3*( g10*cos(theta) + sin(theta)*( g11*cos(phi) + h11*sin(phi) ) ); %(T)
btheta = (Re/r_mag)^3*( g10*sin(theta) - cos(theta)*( g11*cos(phi) + h11*sin(phi) ) ); %(T)
bphi = (Re/r_mag)^3*( g11*sin(phi) - h11*cos(phi) ); %(T)

% Calculate ECI components of magnetic field vector
bi1 = (br*cos(delta) + btheta*sin(delta))*cos(alpha) - bphi*sin(alpha); %(T)
bi2 = (br*cos(delta) + btheta*sin(delta))*sin(alpha) + bphi*cos(alpha); %(T)
bi3 = br*sin(delta) - btheta*cos(delta); %(T)

bi = [bi1 bi2 bi3]'; % Note, this is in the ECI frame!
end



