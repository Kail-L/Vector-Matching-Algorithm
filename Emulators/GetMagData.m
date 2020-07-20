function Ba=GetMagData(ST,pos,t_UTC,t,const)
%GETMAGDATA Gets the Earth Magnetic Field dependent upon the model desired.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Ba] = GETMAGDATA(ST,pos,t,const) takes in position in ECI, along with
% the simulation start time (year,month,day), time since start, and a
% constants structure and determines the magnetic field based on three
% models:
%
% 1.) EarthMagField.m
% 2.) wrldmagm.m
% 3.) igrfmagm.m
%
% The choice of model used is based on how fast the implementation should
% take. The EarthMagField function is a less accurate but speedier
% function, whereas wrldmagm.m is the World Magnetic Model (WMM), used by
% the UK and the U.S. This model is a more accurate yet slower in terms of
% computation. The igrfmagm.m model is the International Geomagnetic
% Reference Field (IGRF Model), but is dependent on MATLAB having the most
% up to date info for the IGRF and may not work. 
%
% SOURCES:
% EarthMagField.m is a modified piece of code created by Professor Ryan 
% Caverly and Richard Forbes, updated by the UMN SmallSat ADCS Team. 
%
% The WMM can be found at:
% https://www.ngdc.noaa.gov/geomag/WMM/DoDWMM.shtml
%
% INPUT PARAMETERS:
% ST = 1x3 array that siginifies the start time of the simulation 
%      ex. ST = [2020 4 20]
% pos = 3x1 position vector in ECI frame, meters
% t_UTC = 1x6 current time in utc format
% t = 1x1 time value of magnetometer sample since sim start, seconds
% const = a constants structure that holds relevant simulation information
%
% OUTPUT PARAMETERS: 
% Ba = 3x1 magnetic field vector in T in the ECI frame.
%
% VARIABLES:
%
% Kail Laughlin
% Updated 6/28/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_a = pos;

if const.MagFieldFlag == 0
    Ba = EarthMagField(r_a,t);
elseif const.MagFieldFlag == 1
    wgs84 = wgs84Ellipsoid('kilometer');
    LLA = eci2lla(r_a',t_UTC);
    
    DCM_EF_E = dcmeci2ecef('IAU-2000/2006',t_UTC);
    DCM_NED_EF = dcmecef2ned(LLA(1),LLA(2));
    DCM_NED_E = DCM_NED_EF*DCM_EF_E;
    
    % Calc mag field in NED frame using WMM %
    [xyz,~,~,~,~] = wrldmagm(LLA(3),LLA(1),LLA(2),decyear(ST)); % nano tesla

    % Convert from NED to ECI %
    Ba = inv(DCM_NED_E)*xyz;
    
    % Convert from nano tesla to tesla
    Ba = Ba*1e-9;

    % make Ba nx3 instead of 3xn
    Ba = Ba';
elseif const.MagFieldFlag == 2
    wgs84 = wgs84Ellipsoid('kilometer');
    LLA = eci2lla(r_a',t_UTC);
    
    DCM_EF_E = dcmeci2ecef('IAU-2000/2006',t_UTC);
    DCM_NED_EF = dcmecef2ned(LLA(1),LLA(2));
    DCM_NED_E = DCM_NED_EF*DCM_EF_E;
    
    % Calc mag field in NED frame using IGRF %   
    [xyz,~,~,~,~,~,~,~,~,~] = igrfmagm(LLA(3),LLA(1),LLA(2),decyear(ST));
    
    % Convert from NED to ECI %
    Ba = inv(DCM_NED_E)*xyz';
    
    % Convert from nano tesla to tesla
    Ba = Ba*1e-9;

    % make Ba nx3 instead of 3xn
    Ba = Ba';
end
end