%% External Data Script %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ExternalData.m processes external data to be used in the VMA.
% This data needs to be in the form of 
%
%       x_out = [r_ECI v_ECI q_S_ECI Om_S_ECI h_S_ECI]
%       ExtVec = [Mag_ECI Mag_S Sun_ECI Sun_S]
%       t = [Gyro/Angular rate time vector]
%       t_sun = [Sun measurement time vector]
%       t_mag = [Mag measurement time vector]
%
%
% where x_out is a n x 16 array of data where n is the length of the data,
% and ExtVec is a nx12 array of measurement data for both sun sensor and
% magnetometer. 
%
% The data is separated into the following variables:
%
% r_ECI - nx3 Inertial position of the spacecraft in m
% v_ECI - nx3 Inertial velocity of the spacecraft in m/s
% q_S_ECI - nx4 True attitude quaternion from ECI frame to S frame
%         - [q0 q1 q2 q3] where q0 is scalar
% Om_S_ECI - nx3 True angular velocity of S frame with respect to inertial ECI
%            frame in rad/s
% h_S_ECI - nx3 True angular momentum vector in kg-m^2/s
% Mag_ECI - nx3 EMF vector modeled in ECI frame in T
% Mag_S - nx3 EMF vector measured in spacecraft frame in T
% Sun_ECI - nx3 sun vector modeled in ECI frame in T
% Sun_S - nx3 sun vector measured in spacecraft frame in T
% t - nx1 gyro time vector in s
%
% Of all this data, the only required components are q_S_ECI and Om_S_ECI.
% To neglect other vectors, replace their values with zeros. The timing
% vectors of t_mag and t_sun are how often the EMF and sun vector are
% measured, respectively. These arrays can be identical to t in the case
% that mag and sun are not used.
%
% Kail Laughlin, 7/7/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SGPB Data %%
% Note: Stanford gravity probe B (SGPB) data is experimental data that
% contains angular rate info of a very high quality gyroscope and
% calibrated magnetometer data. Note that the frames of interest of this
% data are the spacecraft frame S, and the SGPB Inertial frame (SGPBI), NOT
% ECI. 

if const.UseSGPB ==1
    path = pwd;
    addpath([path,'\ExternalData\SGPBData'])

    % Load SGPB data
    load SGPB_Attitude_Data
    load SGPB_Roll_Rate_Data
    load SGPB_Mag_Data
    load SGPB_Mag_Data_Cal
    load SGPB_Sun_Pos_Data

    t = Sens_rr_time_inter;
    t_mag = Mag_time_inter;
    t_sun = GPS_time_inter;

    r_SGPBI = zeros(length(t),3);
    v_SGPBI = zeros(length(t),3);
    psi_S_SGPBI_rg = EA_rg(3,:);
    the_S_SGPBI_rg = EA_rg(2,:);
    phi_S_SGPBI_rg = EA_rg(1,:);
    Euler_angles_rg = [psi_S_SGPBI_rg' the_S_SGPBI_rg' phi_S_SGPBI_rg']; % Store Euler angles
    for lv1=1:length(t)
        q_S_SGPBI(lv1,:) = dcm2quat(eul2dcm(Euler_angles_rg(lv1,:)));  % True quaternion
    end
    Om_S_SGPBI = [rg_true(1,:)' rg_true(2,:)' rg_true(3,:)'];
    h_S_SGPBI = zeros(length(t),3);
    x_out = [r_SGPBI v_SGPBI q_S_SGPBI Om_S_SGPBI h_S_SGPBI];

    B_SGPBI_model = EMF_model_SGPBI';
    B_S = EMF_m_cal';
    Sun_SGPBI = r_SGPB2S_SGPBI_uv';
    Sun_S = r_SGPB2S_SGPB_uv' ;

    ExtVec = [B_SGPBI_model B_S Sun_SGPBI Sun_S];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Alex's HyCUBE Data
% Note: This data is taken from re-entry simulations of HyCUBE ran via Alex
% Hayes' and Ioanis Nompellis' simulation environment. 

if const.UseHyCUBE == 1
    path = pwd;
    addpath([path,'\ExternalData\AlexHyCUBEData'])

    load('simulation_output_5dps_roll.mat')
    t = t_sim';                         % Timespan of sim, seconds
    t_mag = t_sim';                     % Arbitrary time vector for mag needed in post processing
    t_sun = t_sim';                     % Arbitrary time vector for sun needed in post processing
    drl=length(t);                      % Looping variable for integration
    dt = mean(diff(t));                 % Timesteps, seconds
    r_ECI = zeros(length(t),3);
    v_ECI = zeros(length(t),3);

    C_S_ECI = output.DCM_s_ECI;         % DCM of ECI frame to Spacecraft body frame
    q_S_ECI = zeros(length(t),4);       % Quaternion of ECI frame to Spacecraft body frame [scalar vector]
    EulAng_S_ECI = zeros(length(t),3);  % Euler angles of ECI frame to Spacecraft body frame [yaw pitch roll]
    for MC=1:length(C_S_ECI)
        q_S_ECI(MC,:) = dcm2quat(C_S_ECI(:,:,MC));
        if q_S_ECI(MC,1) < 0
            q_S_ECI(MC,:) = -q_S_ECI(MC,:); % Enforce positive scalar
        end
        EulAng_S_ECI(MC,:) = dcm2eul(C_S_ECI(:,:,MC));
    end

    Om_S_ECI = output.w_sECI_s';       % Clean angular rates of spacecraft in body frame, rad/s
    h_S_ECI = zeros(length(t),3);

    x_out = [r_ECI v_ECI q_S_ECI Om_S_ECI h_S_ECI];

    B_ECI_model = zeros(length(t),3);
    B_S = zeros(length(t),3);
    Sun_ECI = zeros(length(t),3);
    Sun_S = zeros(length(t),3);
    ExtVec = [B_ECI_model B_S Sun_ECI Sun_S];
end
