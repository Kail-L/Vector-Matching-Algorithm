%% Earth Mag Field model test
% Note: 
% RefVecs(:,1) = time (s)
% RefVecs(:,2:4) = XYZ coordinates in ECI frame (km)

close all
clear all
clc

load SatPosAtt % load trial data "RefVecs"


XYZ = RefVecs(:,2:4);
time = RefVecs(:,1);

magfield = zeros(length(XYZ),3);
magfield_3rd = zeros(length(XYZ),3);

for ii = 1:length(XYZ)
   
    magfield(ii,:) = EarthMagField(XYZ(ii,:)',time(ii));
    
    magfield_3rd(ii,:) = EarthMagField_3rd_Order(XYZ(ii,:)',time(ii));
    
end

%%
figure(1)
plot3(magfield(:,1),magfield(:,2),magfield(:,3))
grid on
figure(2)
plot3(magfield_3rd(:,1),magfield_3rd(:,2),magfield_3rd(:,3))
grid on
