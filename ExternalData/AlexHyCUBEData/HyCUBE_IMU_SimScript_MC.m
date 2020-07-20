%% HyCUBE IMU Simulation Script %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script takes HyCUBE IMU measurements and propagates an attitude
% solution along with covariance forward in time using only time update
% equations. The goal of the script is to simulate HyCUBE re-entry through
% plasma, a situation where no aiding sensors will be available.
%
%
% Kail Laughlin, Demoz Gebre-Egziabher
% 5/14/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load trajectory data %%
% This data is treated as truth data and is corrupted using GenIMUErr.m

% load('reentry_sim_data.mat')
% load('simulation_output_20hz.mat')
load('simulation_output_5dps_roll.mat')

%% Load covariance data from Kail's Sim %%

load('PinitCon.mat')
load('PinitTac.mat') % actually the initial covariance for NAV grade
load('PinitNav.mat') % actually the initial covariance for TAC grade

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Package Data Into Useable Format %%

MCLength = 30;                      % Number of Monte Carlo runs to perform
d2r = pi/180;                       % degrees to radians conversion
r2d = 1/d2r;                        % radians to degrees conversion
Cbg = eye(3);                       % DCM between spacecraft body and gyro
t = t_sim';                         % Timespan of sim, seconds
drl=length(t);                      % Looping variable for integration
dt = mean(diff(t));                 % Timesteps, seconds
Omega_g_s = output.w_sECI_s';       % Clean angular rates of spacecraft in body frame, rad/s
C_s_ECI = output.DCM_s_ECI;         % DCM of ECI frame to Spacecraft body frame
q_s_ECI = zeros(length(t),4);       % Quaternion of ECI frame to Spacecraft body frame [scalar vector]
EulAng_s_ECI = zeros(length(t),3);  % Euler angles of ECI frame to Spacecraft body frame [yaw pitch roll]
for MC=1:length(C_s_ECI)
    q_s_ECI(MC,:) = dcm2quat(C_s_ECI(:,:,MC));
    if q_s_ECI(MC,1) < 0
        q_s_ECI(MC,:) = -q_s_ECI(MC,:); % Enforce positive scalar
    end
    EulAng_s_ECI(MC,:) = dcm2eul(C_s_ECI(:,:,MC));
end

% Place angular rate data into format for GenIMUerr
imu_good=[t Omega_g_s(:,1) Omega_g_s(:,2) Omega_g_s(:,3)...
    zeros(length(Omega_g_s),3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
direc = pwd;
font_size = 15;
line_size = 15;
line_width = 1;
TruCol = 'b';   % Color of truth data
ConCol = 'r';   % Color of CON data
TacCol = 'g';   % Color of TAC data
NavCol = 'k';   % Color of NAV data
SaveLoc = [direc,'\MCFigures']; % Location (filepath) to save figures to
SM = 1;                       % Number of sigma to plot (1-sigma, 2-sigma, etc.)
plotTime = 1;                 % Boolean to plot in seconds or hours (0 - seconds, 1 - hours)
if plotTime==2          
    t_plot = t./3600;
    timeval='hrs';
elseif plotTime==1
    t_plot=t./60;
    timeval='min';
elseif plotTije==0
    t_plot = t;
    timeval='s';
end
XLimit = [0 t_plot(end)];          % Xlimit on plotting
YLimitEAngCon = [-30 30];             % Ylimit on euler angle error plot
YLimitEAngTacNav = [-5 5];             % Ylimit on euler angle error plot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MC Storage arrays %%
% Noisy Quaternions %
q_s_hat_ECI0_Con_MC = zeros(drl,MCLength);
q_s_hat_ECI1_Con_MC = zeros(drl,MCLength);
q_s_hat_ECI2_Con_MC = zeros(drl,MCLength);
q_s_hat_ECI3_Con_MC = zeros(drl,MCLength);
q_s_hat_ECI0_Tac_MC = zeros(drl,MCLength);
q_s_hat_ECI1_Tac_MC = zeros(drl,MCLength);
q_s_hat_ECI2_Tac_MC = zeros(drl,MCLength);
q_s_hat_ECI3_Tac_MC = zeros(drl,MCLength);
q_s_hat_ECI0_Nav_MC = zeros(drl,MCLength);
q_s_hat_ECI1_Nav_MC = zeros(drl,MCLength);
q_s_hat_ECI2_Nav_MC = zeros(drl,MCLength);
q_s_hat_ECI3_Nav_MC = zeros(drl,MCLength);

% Error Quaternions %
q_s_s_hat0_Con_MC = zeros(drl,MCLength);
q_s_s_hat1_Con_MC = zeros(drl,MCLength);
q_s_s_hat2_Con_MC = zeros(drl,MCLength);
q_s_s_hat3_Con_MC = zeros(drl,MCLength);    
q_s_s_hat0_Tac_MC = zeros(drl,MCLength);
q_s_s_hat1_Tac_MC = zeros(drl,MCLength);
q_s_s_hat2_Tac_MC = zeros(drl,MCLength);
q_s_s_hat3_Tac_MC = zeros(drl,MCLength);
q_s_s_hat0_Nav_MC = zeros(drl,MCLength);
q_s_s_hat1_Nav_MC = zeros(drl,MCLength);
q_s_s_hat2_Nav_MC = zeros(drl,MCLength);
q_s_s_hat3_Nav_MC = zeros(drl,MCLength);


% "Noisy" Euler Angles %
psi_s_ECI_Con_MC = zeros(drl,MCLength);
the_s_ECI_Con_MC = zeros(drl,MCLength);
phi_s_ECI_Con_MC = zeros(drl,MCLength);
psi_s_ECI_Nav_MC = zeros(drl,MCLength);
the_s_ECI_Nav_MC = zeros(drl,MCLength);
phi_s_ECI_Nav_MC = zeros(drl,MCLength);
psi_s_ECI_Tac_MC = zeros(drl,MCLength);
the_s_ECI_Tac_MC = zeros(drl,MCLength);
phi_s_ECI_Tac_MC = zeros(drl,MCLength);

% Error Euler Angles %
psi_Err_Con_MC = zeros(drl,MCLength);
the_Err_Con_MC = zeros(drl,MCLength);
phi_Err_Con_MC = zeros(drl,MCLength);
psi_Err_Nav_MC = zeros(drl,MCLength);
the_Err_Nav_MC = zeros(drl,MCLength);
phi_Err_Nav_MC = zeros(drl,MCLength);
psi_Err_Tac_MC = zeros(drl,MCLength);
the_Err_Tac_MC = zeros(drl,MCLength);
phi_Err_Tac_MC = zeros(drl,MCLength);

% Error between 3-axis of s and s_hat frames %
epsilon_Con_MC = zeros(drl,MCLength);
epsilon_Tac_MC = zeros(drl,MCLength);
epsilon_Nav_MC = zeros(drl,MCLength);

% Standard deviations of states %
std_qe1_Con_MC = zeros(drl,MCLength);
std_qe1_Tac_MC = zeros(drl,MCLength);
std_qe1_Nav_MC = zeros(drl,MCLength);
std_qe2_Con_MC = zeros(drl,MCLength);
std_qe2_Tac_MC = zeros(drl,MCLength);
std_qe2_Nav_MC = zeros(drl,MCLength);
std_qe3_Con_MC = zeros(drl,MCLength);
std_qe3_Tac_MC = zeros(drl,MCLength);
std_qe3_Nav_MC = zeros(drl,MCLength);
std_psi_e_Con_MC = zeros(drl,MCLength);
std_psi_e_Tac_MC = zeros(drl,MCLength);
std_psi_e_Nav_MC = zeros(drl,MCLength);
std_the_e_Con_MC = zeros(drl,MCLength);
std_the_e_Tac_MC = zeros(drl,MCLength);
std_the_e_Nav_MC = zeros(drl,MCLength);
std_phi_e_Con_MC = zeros(drl,MCLength);
std_phi_e_Tac_MC = zeros(drl,MCLength);
std_phi_e_Nav_MC = zeros(drl,MCLength);
std_epsilon_Con_MC = zeros(drl,MCLength);
std_epsilon_Tac_MC = zeros(drl,MCLength);
std_epsilon_Nav_MC = zeros(drl,MCLength);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Monte Carlo %%
for MC = 1:MCLength
    %% Corrupt IMU Data Dependent on Grade of IMU %%
    disp(['MC Run: ',num2str(MC),'/',num2str(MCLength)])

    % Conventional IMU Grade %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IRb_Con = 8.7266e-4;
    tau_Con = 300;
    OutNoise_Con = 8.7266e-4;
    [imu_corrupt_Con,~] = emulateGyroHyCUBE(t,imu_good(:,2:4),IRb_Con,tau_Con,OutNoise_Con);

    % Re-Package Data %
    IMU_Conventional=[imu_good(:,1:4) imu_corrupt_Con];

    % Tactical IMU Grade %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IRb_Tac = 1.6968e-6;
    tau_Tac  = 1000;
    OutNoise_Tac  = 2.9671e-5;
    [imu_corrupt_Tac ,~] = emulateGyroHyCUBE(t,imu_good(:,2:4),IRb_Tac ,tau_Tac ,OutNoise_Tac );

    % Re-Package Data %
    IMU_Tactical=[imu_good(:,1:4) imu_corrupt_Tac];

    % Navigation IMU Grade %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    IRb_Nav = 1.4544e-8;
    tau_Nav  = 3600;
    OutNoise_Nav  = 1.3963e-5;
    [imu_corrupt_Nav ,~] = emulateGyroHyCUBE(t,imu_good(:,2:4),IRb_Nav ,tau_Nav ,OutNoise_Nav );

    % Re-Package Data %
    IMU_Navigation=[imu_good(:,1:4) imu_corrupt_Nav];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% EKF Time Update Only %%
    % Set up Data and initialize matrices for  s p e e d %

    % "Clean" IMU measurements
    w_s_Con_clean = IMU_Conventional(:,2:4);
    w_s_Nav_clean = IMU_Navigation(:,2:4);
    w_s_Tac_clean = IMU_Tactical(:,2:4);
    % "Noisy" IMU measurements
    w_s_Con_noisy = IMU_Conventional(:,5:7);
    w_s_Nav_noisy = IMU_Navigation(:,5:7);
    w_s_Tac_noisy = IMU_Tactical(:,5:7);
    % Difference between true and noisy IMU measurements %
    w_s_Con_diff = w_s_Con_noisy-w_s_Con_clean;
    w_s_Nav_diff = w_s_Nav_noisy-w_s_Nav_clean;
    w_s_Tac_diff = w_s_Tac_noisy-w_s_Tac_clean;
    % Noisy Quaternions %
    q_s_hat_ECI_Con = zeros(drl,4);
    q_s_hat_ECI_Tac = zeros(drl,4);
    q_s_hat_ECI_Nav = zeros(drl,4);
    % Error Quaternions %
    q_s_s_hat_Con = zeros(drl,4);
    q_s_s_hat_Tac = zeros(drl,4);
    q_s_s_hat_Nav = zeros(drl,4);
    % True Euler Angles %
    psi_s_ECI = EulAng_s_ECI(:,1);
    the_s_ECI = EulAng_s_ECI(:,2);
    phi_s_ECI = EulAng_s_ECI(:,3);
    % "Noisy" Euler Angles %
    psi_s_ECI_Con = zeros(drl,1);
    the_s_ECI_Con = zeros(drl,1);
    phi_s_ECI_Con = zeros(drl,1);
    psi_s_ECI_Nav = zeros(drl,1);
    the_s_ECI_Nav = zeros(drl,1);
    phi_s_ECI_Nav = zeros(drl,1);
    psi_s_ECI_Tac = zeros(drl,1);
    the_s_ECI_Tac = zeros(drl,1);
    phi_s_ECI_Tac = zeros(drl,1);
    % Error Euler Angles %
    psi_Err_Con = zeros(drl,1);
    the_Err_Con = zeros(drl,1);
    phi_Err_Con = zeros(drl,1);
    psi_Err_Nav = zeros(drl,1);
    the_Err_Nav = zeros(drl,1);
    phi_Err_Nav = zeros(drl,1);
    psi_Err_Tac = zeros(drl,1);
    the_Err_Tac = zeros(drl,1);
    phi_Err_Tac = zeros(drl,1);
    % Error between 3-axis of s and s_hat frames %
    epsilon_Con = zeros(drl,1);
    epsilon_Tac = zeros(drl,1);
    epsilon_Nav = zeros(drl,1);
    % State Covariance %
    P_stor_Con = zeros(drl,6);
    P_stor_Tac = zeros(drl,6);
    P_stor_Nav = zeros(drl,6);
    P_dq_Con = zeros(4,4);     % Covariance matrix used to convert quat covariance to eul covariance
    P_dq_Tac = zeros(4,4);     % Covariance matrix used to convert quat covariance to eul covariance
    P_dq_Nav = zeros(4,4);     % Covariance matrix used to convert quat covariance to eul covariance
    % Standard deviations of states %
    std_qe1_Con = zeros(drl,1);
    std_qe1_Tac = zeros(drl,1);
    std_qe1_Nav = zeros(drl,1);
    std_qe2_Con = zeros(drl,1);
    std_qe2_Tac = zeros(drl,1);
    std_qe2_Nav = zeros(drl,1);
    std_qe3_Con = zeros(drl,1);
    std_qe3_Tac = zeros(drl,1);
    std_qe3_Nav = zeros(drl,1);
    std_psi_e_Con = zeros(drl,1);
    std_psi_e_Tac = zeros(drl,1);
    std_psi_e_Nav = zeros(drl,1);
    std_the_e_Con = zeros(drl,1);
    std_the_e_Tac = zeros(drl,1);
    std_the_e_Nav = zeros(drl,1);
    std_phi_e_Con = zeros(drl,1);
    std_phi_e_Tac = zeros(drl,1);
    std_phi_e_Nav = zeros(drl,1);
    std_epsilon_Con = zeros(drl,1);
    std_epsilon_Tac = zeros(drl,1);
    std_epsilon_Nav = zeros(drl,1);
    
    clear 'IMUDat_Conventional.mat'
    clear 'IMUDat_Tactical.mat'
    clear 'IMUDat_Navigation.mat'
    %% EKF Initialization %%

    % Initial Covariance Values %
    % PinitCon,PinitTac, and PinitNav taken as state covariance after Kail's Sim converges
    %
    % PinitCon created from running with Con IMU characteristics
    % taken from GenIMUErr.

    P0_Con=PinitCon;
    P0_Tac=PinitNav; % Note, the values in PinitNav actually correspond to the tactical IMU
    P0_Nav=PinitTac; % Note, the values in the PinitTac actually correspond to the navigational IMU

    P_stor_Con(1,:) = [P0_Con(1,1) P0_Con(2,2) P0_Con(3,3) P0_Con(4,4) P0_Con(5,5) P0_Con(6,6)];        
    P_stor_Tac(1,:) = [P0_Tac(1,1) P0_Tac(2,2) P0_Tac(3,3) P0_Tac(4,4) P0_Tac(5,5) P0_Tac(6,6)];  
    P_stor_Nav(1,:) = [P0_Nav(1,1) P0_Nav(2,2) P0_Nav(3,3) P0_Nav(4,4) P0_Nav(5,5) P0_Nav(6,6)];  

    % Initial Standard Deviations on Error States %
    std_qe1_Con(1) = sqrt(P0_Con(1,1));
    std_qe1_Tac(1) = sqrt(P0_Tac(1,1));
    std_qe1_Nav(1) = sqrt(P0_Nav(1,1));
    std_qe2_Con(1) = sqrt(P0_Con(2,2));
    std_qe2_Tac(1) = sqrt(P0_Tac(2,2));
    std_qe2_Nav(1) = sqrt(P0_Nav(2,2));
    std_qe3_Con(1) = sqrt(P0_Con(3,3));
    std_qe3_Tac(1) = sqrt(P0_Tac(3,3));
    std_qe3_Nav(1) = sqrt(P0_Nav(3,3));

    % Initial Process Noise Covariance Matrices %
    Qw_Con = [(OutNoise_Con)^2*eye(3) zeros(3,3);     % Gyro Output noise
            zeros(3,3) 2*(IRb_Con)^2/tau_Con*eye(3)]; % Gyro in-run bias repeat.

    Qw_Tac = [(OutNoise_Tac)^2*eye(3) zeros(3,3);     % Gyro Output noise
            zeros(3,3) 2*(IRb_Tac)^2/tau_Tac*eye(3)]; % Gyro in-run bias repeat.

    Qw_Nav = [(OutNoise_Nav)^2*eye(3) zeros(3,3);     % Gyro Output noise
            zeros(3,3) 2*(IRb_Nav)^2/tau_Nav*eye(3)]; % Gyro in-run bias repeat.

    % Initial Quaternion Estimate (Use true quaterion) %
    q_s_hat_ECI_Con(1,:) = q_s_ECI(1,:);
    q_s_hat_ECI_Tac(1,:) = q_s_ECI(1,:);
    q_s_hat_ECI_Nav(1,:) = q_s_ECI(1,:);

    % Determine Initial Angular Error Between S3 and S3_hat axes %%%%%%%%%%%%%%

    % Estimate DCM from ECI to s frame
    CsECI_hat_0_Con = quat2dcm(q_s_hat_ECI_Con(1,:));
    CsECI_hat_0_Tac = quat2dcm(q_s_hat_ECI_Tac(1,:));
    CsECI_hat_0_Nav = quat2dcm(q_s_hat_ECI_Nav(1,:));

    % Error DCM between s and s_hat frame
    Css_hat_0_Con = quat2dcm(q_s_ECI(1,:))*inv(CsECI_hat_0_Con);
    Css_hat_0_Tac = quat2dcm(q_s_ECI(1,:))*inv(CsECI_hat_0_Tac);
    Css_hat_0_Nav = quat2dcm(q_s_ECI(1,:))*inv(CsECI_hat_0_Nav);

    % Initial Euler Angles (rad) %
    E_s_ECI_hat_0_Con = dcm2eul(quat2dcm(q_s_hat_ECI_Con(1,:)));
    E_s_ECI_hat_0_Tac = dcm2eul(quat2dcm(q_s_hat_ECI_Tac(1,:)));
    E_s_ECI_hat_0_Nav = dcm2eul(quat2dcm(q_s_hat_ECI_Nav(1,:)));

    psi_s_ECI_Con(1) = E_s_ECI_hat_0_Con(1);
    the_s_ECI_Con(1) = E_s_ECI_hat_0_Con(2);
    phi_s_ECI_Con(1) = E_s_ECI_hat_0_Con(3);
    psi_s_ECI_Tac(1) = E_s_ECI_hat_0_Tac(1);
    the_s_ECI_Tac(1) = E_s_ECI_hat_0_Tac(2);
    phi_s_ECI_Tac(1) = E_s_ECI_hat_0_Tac(3);
    psi_s_ECI_Nav(1) = E_s_ECI_hat_0_Nav(1);
    the_s_ECI_Nav(1) = E_s_ECI_hat_0_Nav(2);
    phi_s_ECI_Nav(1) = E_s_ECI_hat_0_Nav(3);

    vE = [0 0 1]';               % 3-axis of s frame
    vS_Con = Css_hat_0_Con*vE;   % 3-axis of s_hat frame
    vS_Tac = Css_hat_0_Tac*vE;   % 3-axis of s_hat frame
    vS_Nav = Css_hat_0_Nav*vE;   % 3-axis of s_hat frame

    % Initial Angular error between 3-axis of s and s_hat frames %
    epsilon_Con(1) = acos(dot(vS_Con,vE)/(norm(vS_Con)*norm(vE)))*180/pi;
    epsilon_Tac(1) = acos(dot(vS_Tac,vE)/(norm(vS_Tac)*norm(vE)))*180/pi;
    epsilon_Nav(1) = acos(dot(vS_Nav,vE)/(norm(vS_Nav)*norm(vE)))*180/pi;

    % Initial Quaternion Error between s and s_hat %
    q_s_s_hat_Con(1,:) = quatmult(q_s_ECI(1,:)',quatInverse(q_s_hat_ECI_Con(1,:)));
    q_s_s_hat_Tac(1,:) = quatmult(q_s_ECI(1,:)',quatInverse(q_s_hat_ECI_Tac(1,:)));
    q_s_s_hat_Nav(1,:) = quatmult(q_s_ECI(1,:)',quatInverse(q_s_hat_ECI_Nav(1,:)));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Run EKF (Time Update Only) %%
    for k = 2:drl
        w_hat_Con = w_s_Con_noisy(k,:)';  % Estimate angular velocity, body frame (no bias subracted)
        F_Con = [-sk(w_hat_Con) -1/2*eye(3); 
                 zeros(3) (-1/tau_Con)*eye(3)];
        G_Con = [-1/2*Cbg zeros(3);         
                 zeros(3) eye(3)];

        w_hat_Tac = w_s_Tac_noisy(k,:)';  % Estimate angular velocity, body frame (no bias subracted)
        F_Tac = [-sk(w_hat_Tac) -1/2*eye(3); 
                 zeros(3) (-1/tau_Tac)*eye(3)];
        G_Tac = [-1/2*Cbg zeros(3);         
                 zeros(3) eye(3)];

        w_hat_Nav = w_s_Nav_noisy(k,:)';  % Estimate angular velocity, body frame (no bias subracted)
        F_Nav = [-sk(w_hat_Nav) -1/2*eye(3); 
                 zeros(3) (-1/tau_Nav)*eye(3)];
        G_Nav = [-1/2*Cbg zeros(3);         
                 zeros(3) eye(3)];         

        % Covariance Update %
        Phi_Con = expm(F_Con*dt);
        Phi_Tac = expm(F_Tac*dt);
        Phi_Nav = expm(F_Nav*dt);

        Cd_Con = disrw(F_Con,G_Con,dt,Qw_Con);
        Cd_Tac = disrw(F_Tac,G_Tac,dt,Qw_Tac);    
        Cd_Nav = disrw(F_Nav,G_Nav,dt,Qw_Nav);  

        P_Con = Phi_Con*P0_Con*Phi_Con' + Cd_Con;
        P_Tac = Phi_Tac*P0_Tac*Phi_Tac' + Cd_Tac;
        P_Nav = Phi_Nav*P0_Nav*Phi_Nav' + Cd_Nav;

        % Store P %
        P_stor_Con(k,:) = [P_Con(1,1) P_Con(2,2) P_Con(3,3) P_Con(4,4) P_Con(5,5) P_Con(6,6)];        
        P_stor_Tac(k,:) = [P_Tac(1,1) P_Tac(2,2) P_Tac(3,3) P_Tac(4,4) P_Tac(5,5) P_Tac(6,6)];  
        P_stor_Nav(k,:) = [P_Nav(1,1) P_Nav(2,2) P_Nav(3,3) P_Nav(4,4) P_Nav(5,5) P_Nav(6,6)]; 
        P0_Con = P_Con;
        P0_Tac = P_Tac;
        P0_Nav = P_Nav;

        % Estimate Quaternion Update %
        q_s_hat_ECI_Con(k,:) = quatIntegrate((w_hat_Con)',q_s_hat_ECI_Con(k-1,:)',dt);
        q_s_hat_ECI_Tac(k,:) = quatIntegrate((w_hat_Tac)',q_s_hat_ECI_Tac(k-1,:)',dt);
        q_s_hat_ECI_Nav(k,:) = quatIntegrate((w_hat_Nav)',q_s_hat_ECI_Nav(k-1,:)',dt);

        % Enforce Positive Scalar on Quaternion %
        if q_s_hat_ECI_Con(k,1) < 0
            q_s_hat_ECI_Con(k,:) = -q_s_hat_ECI_Con(k,:);
        end

        if q_s_hat_ECI_Tac(k,1) < 0
            q_s_hat_ECI_Tac(k,:) = -q_s_hat_ECI_Tac(k,:);
        end    

        if q_s_hat_ECI_Nav(k,1) < 0
            q_s_hat_ECI_Nav(k,:) = -q_s_hat_ECI_Nav(k,:);
        end

        % Updated DCM %
        C_s_hat_ECI_Con = quat2dcm(q_s_hat_ECI_Con(k,:));
        C_s_hat_ECI_Tac = quat2dcm(q_s_hat_ECI_Tac(k,:));
        C_s_hat_ECI_Nav = quat2dcm(q_s_hat_ECI_Nav(k,:));

        % Error Quaternion Update %
        C_s_s_hat_Con = C_s_ECI(:,:,k)*(C_s_hat_ECI_Con)';
        C_s_s_hat_Tac = C_s_ECI(:,:,k)*(C_s_hat_ECI_Tac)';
        C_s_s_hat_Nav = C_s_ECI(:,:,k)*(C_s_hat_ECI_Nav)';

        q_s_s_hat_Con(k,:) = dcm2quat(C_s_s_hat_Con);
        q_s_s_hat_Tac(k,:) = dcm2quat(C_s_s_hat_Tac);
        q_s_s_hat_Nav(k,:) = dcm2quat(C_s_s_hat_Nav);

        % Angular Error %
        vS_Con = C_s_s_hat_Con*vE;
        vS_Tac = C_s_s_hat_Tac*vE;
        vS_Nav = C_s_s_hat_Nav*vE;

        epsilon_Con(k) = acos(dot(vS_Con,vE)/(norm(vS_Con)*norm(vE)))*180/pi;
        epsilon_Tac(k) = acos(dot(vS_Tac,vE)/(norm(vS_Tac)*norm(vE)))*180/pi;
        epsilon_Nav(k) = acos(dot(vS_Nav,vE)/(norm(vS_Nav)*norm(vE)))*180/pi;

        % STD Deviation of error quaternions %
        std_qe1_Con(k) = sqrt(P_Con(1,1));
        std_qe2_Con(k) = sqrt(P_Con(2,2));
        std_qe3_Con(k) = sqrt(P_Con(3,3));
        std_qe1_Tac(k) = sqrt(P_Tac(1,1));
        std_qe2_Tac(k) = sqrt(P_Tac(2,2));
        std_qe3_Tac(k) = sqrt(P_Tac(3,3));    
        std_qe1_Nav(k) = sqrt(P_Nav(1,1));
        std_qe2_Nav(k) = sqrt(P_Nav(2,2));
        std_qe3_Nav(k) = sqrt(P_Nav(3,3)); 

        % Estimated Euler Angles %
        Eul_s_hat_ECI_Con = quat2eul(q_s_hat_ECI_Con(k,:));
        Eul_s_hat_ECI_Tac = quat2eul(q_s_hat_ECI_Tac(k,:));
        Eul_s_hat_ECI_Nav = quat2eul(q_s_hat_ECI_Nav(k,:));

        psi_s_ECI_Con(k) = Eul_s_hat_ECI_Con(1);
        the_s_ECI_Con(k) = Eul_s_hat_ECI_Con(2);
        phi_s_ECI_Con(k) = Eul_s_hat_ECI_Con(3);    
        psi_s_ECI_Tac(k) = Eul_s_hat_ECI_Tac(1);
        the_s_ECI_Tac(k) = Eul_s_hat_ECI_Tac(2);
        phi_s_ECI_Tac(k) = Eul_s_hat_ECI_Tac(3);  
        psi_s_ECI_Nav(k) = Eul_s_hat_ECI_Nav(1);
        the_s_ECI_Nav(k) = Eul_s_hat_ECI_Nav(2);
        phi_s_ECI_Nav(k) = Eul_s_hat_ECI_Nav(3); 

        % Estimated Euler Angle Errors %
        Eul_Err_Con = dcm2eul(C_s_s_hat_Con);
        Eul_Err_Tac = dcm2eul(C_s_s_hat_Tac);
        Eul_Err_Nav = dcm2eul(C_s_s_hat_Nav);


        psi_Err_Con(k) = Eul_Err_Con(1);
        the_Err_Con(k) = Eul_Err_Con(2);
        phi_Err_Con(k) = Eul_Err_Con(3);    
        psi_Err_Tac(k) = Eul_Err_Tac(1);
        the_Err_Tac(k) = Eul_Err_Tac(2);
        phi_Err_Tac(k) = Eul_Err_Tac(3);  
        psi_Err_Nav(k) = Eul_Err_Nav(1);
        the_Err_Nav(k) = Eul_Err_Nav(2);
        phi_Err_Nav(k) = Eul_Err_Nav(3);     

        % STD of Euler Angle Errors %

        % Note to Demoz: Take a look at the QuatJoc.m function. I am pretty
        % positive it is correct, but there may be some finnicky details due to
        % the fact that we only have a covariance on the vector portion of the
        % quaternion.

        Gq_Con = QuatJoc(q_s_s_hat_Con(k,:)');
        Gq_Tac = QuatJoc(q_s_s_hat_Tac(k,:)');
        Gq_Nav = QuatJoc(q_s_s_hat_Nav(k,:)');

        % Note: The above jacobians are with respect to a full quaternion (4x1)
        % whereas our error quaternion is only 3x1, as we assume the scalar is
        % constant. As such, we need to augment either the covariance matrix to
        % include a covariance of 0 on the scalar portion, or remove the first
        % column of the jacobian to remove the calculation on the scalar
        % component. For now, I choose to augment the Jacobian.

        P_dq_Con = P_Con(1:3,1:3);
        P_dq_Tac = P_Tac(1:3,1:3);
        P_dq_Nav = P_Nav(1:3,1:3);
        Peul_Con = Gq_Con(:,2:4)*P_dq_Con*Gq_Con(:,2:4)';
        Peul_Tac = Gq_Tac(:,2:4)*P_dq_Tac*Gq_Tac(:,2:4)';
        Peul_Nav = Gq_Nav(:,2:4)*P_dq_Nav*Gq_Nav(:,2:4)';
        std_psi_e_Con(k) = sqrt(Peul_Con(1,1));
        std_psi_e_Tac(k) = sqrt(Peul_Tac(1,1));
        std_psi_e_Nav(k) = sqrt(Peul_Nav(1,1));
        std_the_e_Con(k) = sqrt(Peul_Con(2,2));
        std_the_e_Tac(k) = sqrt(Peul_Tac(2,2));
        std_the_e_Nav(k) = sqrt(Peul_Nav(2,2));
        std_phi_e_Con(k) = sqrt(Peul_Con(3,3));
        std_phi_e_Tac(k) = sqrt(Peul_Tac(3,3));
        std_phi_e_Nav(k) = sqrt(Peul_Nav(3,3));

        % STD of Epsilon %
        % Use a similar Jacobian method used in Euler angle STD %
        Eq_Con = EpsJoc(q_s_s_hat_Con);
        Eq_Tac = EpsJoc(q_s_s_hat_Tac);
        Eq_Nav = EpsJoc(q_s_s_hat_Nav);
        P_eps_Con=P_Con(1:3,1:3);
        P_eps_Tac=P_Tac(1:3,1:3);
        P_eps_Nav=P_Nav(1:3,1:3);
        Peps_Con=Eq_Con(2:4)*P_eps_Con*Eq_Con(2:4)';
        Peps_Tac=Eq_Tac(2:4)*P_eps_Tac*Eq_Tac(2:4)';
        Peps_Nav=Eq_Nav(2:4)*P_eps_Nav*Eq_Nav(2:4)';
        std_epsilon_Con(k) = sqrt(Peps_Con)*r2d;
        std_epsilon_Tac(k) = sqrt(Peps_Tac)*r2d;
        std_epsilon_Nav(k) = sqrt(Peps_Nav)*r2d;

    end
    %% MC Storage
    % Noisy Quaternions %
    q_s_hat_ECI0_Con_MC(:,MC) = q_s_hat_ECI_Con(:,1);
    q_s_hat_ECI1_Con_MC(:,MC)  = q_s_hat_ECI_Con(:,2);
    q_s_hat_ECI2_Con_MC(:,MC)  = q_s_hat_ECI_Con(:,3);
    q_s_hat_ECI3_Con_MC(:,MC)  = q_s_hat_ECI_Con(:,4);
    q_s_hat_ECI0_Tac_MC(:,MC)  = q_s_hat_ECI_Tac(:,1);
    q_s_hat_ECI1_Tac_MC(:,MC)  = q_s_hat_ECI_Tac(:,2);
    q_s_hat_ECI2_Tac_MC(:,MC)  = q_s_hat_ECI_Tac(:,3);
    q_s_hat_ECI3_Tac_MC(:,MC)  = q_s_hat_ECI_Tac(:,4);
    q_s_hat_ECI0_Nav_MC(:,MC)  = q_s_hat_ECI_Nav(:,1);
    q_s_hat_ECI1_Nav_MC(:,MC)  = q_s_hat_ECI_Nav(:,2);
    q_s_hat_ECI2_Nav_MC(:,MC)  = q_s_hat_ECI_Nav(:,3);
    q_s_hat_ECI3_Nav_MC(:,MC)  = q_s_hat_ECI_Nav(:,4);

    % Error Quaternions %
    q_s_s_hat0_Con_MC(:,MC)  = q_s_s_hat_Con(:,1);
    q_s_s_hat1_Con_MC(:,MC)  = q_s_s_hat_Con(:,2);
    q_s_s_hat2_Con_MC(:,MC)  = q_s_s_hat_Con(:,3);
    q_s_s_hat3_Con_MC(:,MC)  = q_s_s_hat_Con(:,4);   
    q_s_s_hat0_Tac_MC(:,MC)  = q_s_s_hat_Tac(:,1);
    q_s_s_hat1_Tac_MC(:,MC)  = q_s_s_hat_Tac(:,2);
    q_s_s_hat2_Tac_MC(:,MC)  = q_s_s_hat_Tac(:,3);
    q_s_s_hat3_Tac_MC(:,MC)  = q_s_s_hat_Tac(:,4);
    q_s_s_hat0_Nav_MC(:,MC)  = q_s_s_hat_Nav(:,1);
    q_s_s_hat1_Nav_MC(:,MC)  = q_s_s_hat_Nav(:,2);
    q_s_s_hat2_Nav_MC(:,MC)  = q_s_s_hat_Nav(:,3);
    q_s_s_hat3_Nav_MC(:,MC)  = q_s_s_hat_Nav(:,4);


    % "Noisy" Euler Angles %
    psi_s_ECI_Con_MC(:,MC)  = psi_s_ECI_Con;
    the_s_ECI_Con_MC(:,MC)  = the_s_ECI_Con;
    phi_s_ECI_Con_MC(:,MC)  = phi_s_ECI_Con;
    psi_s_ECI_Nav_MC(:,MC)  = psi_s_ECI_Nav;
    the_s_ECI_Nav_MC(:,MC)  = the_s_ECI_Nav;
    phi_s_ECI_Nav_MC(:,MC)  = phi_s_ECI_Nav;
    psi_s_ECI_Tac_MC(:,MC)  = psi_s_ECI_Tac;
    the_s_ECI_Tac_MC(:,MC)  = the_s_ECI_Tac;
    phi_s_ECI_Tac_MC(:,MC)  = phi_s_ECI_Tac;

    % Error Euler Angles %
    psi_Err_Con_MC(:,MC)  = psi_Err_Con;
    the_Err_Con_MC(:,MC)  = the_Err_Con;
    phi_Err_Con_MC(:,MC)  = phi_Err_Con;
    psi_Err_Nav_MC(:,MC)  = psi_Err_Nav;
    the_Err_Nav_MC(:,MC)  = the_Err_Nav;
    phi_Err_Nav_MC(:,MC)  = phi_Err_Nav;
    psi_Err_Tac_MC(:,MC)  = psi_Err_Tac;
    the_Err_Tac_MC(:,MC)  = the_Err_Tac;
    phi_Err_Tac_MC(:,MC)  = phi_Err_Tac;

    % Error between 3-axis of s and s_hat frames %
    epsilon_Con_MC(:,MC)  = epsilon_Con;
    epsilon_Tac_MC(:,MC)  = epsilon_Tac;
    epsilon_Nav_MC(:,MC)  = epsilon_Nav;

    % Standard deviations of states %
    std_qe1_Con_MC(:,MC)  = std_qe1_Con;
    std_qe1_Tac_MC(:,MC)  = std_qe1_Tac;
    std_qe1_Nav_MC(:,MC)  = std_qe1_Nav;
    std_qe2_Con_MC(:,MC)  = std_qe2_Con;
    std_qe2_Tac_MC(:,MC)  = std_qe2_Tac;
    std_qe2_Nav_MC(:,MC)  = std_qe2_Nav;
    std_qe3_Con_MC(:,MC)  = std_qe3_Con;
    std_qe3_Tac_MC(:,MC)  = std_qe3_Tac;
    std_qe3_Nav_MC(:,MC)  = std_qe3_Nav;
    std_psi_e_Con_MC(:,MC)  = std_psi_e_Con;
    std_psi_e_Tac_MC(:,MC)  = std_psi_e_Tac;
    std_psi_e_Nav_MC(:,MC)  = std_psi_e_Nav;
    std_the_e_Con_MC(:,MC)  = std_the_e_Con;
    std_the_e_Tac_MC(:,MC)  = std_the_e_Tac;
    std_the_e_Nav_MC(:,MC)  = std_the_e_Nav;
    std_phi_e_Con_MC(:,MC)  = std_phi_e_Con;
    std_phi_e_Tac_MC(:,MC)  = std_phi_e_Tac;
    std_phi_e_Nav_MC(:,MC)  = std_phi_e_Nav;
    std_epsilon_Con_MC(:,MC) = std_epsilon_Con;
    std_epsilon_Tac_MC(:,MC) = std_epsilon_Tac;
    std_epsilon_Nav_MC(:,MC) = std_epsilon_Nav;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Results %%
% Font size, line size, and line width %
font_size = 18;
yfont_size = 20;
leg_font_size = 12;
axes_font_size = 12;
line_size = 15;
line_width = 1;

% Angular Rates of NAV grade %
if MCLength==1
    Col = 'r';
    axFont = 16;
%     subplot(321);
%     plot(t,r2d*imu_good(:,2),'b');
%     grid;ylabel('$$\omega^{SE}_{S1}$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
%     subplot(322);
%     plot(t,r2d*imu_corrupt_Nav(:,1),Col);
%     grid;ylabel('$$\omega^{SE}_{S1,m}$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
%     subplot(323);
%     plot(t,r2d*imu_good(:,3),'b');
%     grid;ylabel('$$\omega^{SE}_{S2}$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
%     subplot(324);
%     plot(t,r2d*imu_corrupt_Nav(:,2),Col);
%     grid;ylabel('$$\omega^{SE}_{S2,m}$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
%     subplot(325);
%     plot(t,r2d*imu_good(:,4),'b');
%     grid;ylabel('$$\omega^{SE}_{S3}$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
%     subplot(326)
%     plot(t,r2d*imu_corrupt_Nav(:,3),Col);
%     grid;ylabel('$$\omega^{SE}_{S2,m}$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
    subplot(311);
    plot(t,r2d*imu_corrupt_Nav(:,1),Col);
    grid;ylabel('$$\omega^{SE}_{S1,m}$$','Interpreter','Latex','FontSize',yfont_size);
    ax = gca;
    ax.FontSize = axFont;
    title('$$\omega^{SE}$$ (deg/s)','Interpreter','Latex','FontSize',yfont_size);
    subplot(312);
    plot(t,r2d*imu_corrupt_Nav(:,2),Col);
    grid;ylabel('$$\omega^{SE}_{S2,m}$$','Interpreter','Latex','FontSize',yfont_size);
    ax = gca;
    ax.FontSize = axFont;
    subplot(313)
    plot(t,r2d*imu_corrupt_Nav(:,3),Col);
    grid;ylabel('$$\omega^{SE}_{S2,m}$$','Interpreter','Latex','FontSize',yfont_size);
    ax = gca;
    ax.FontSize = axFont;
end

% CON Euler Angle Errors %
fig=gcf;
figure(fig.Number+1)
EulErrAx(1) = subplot(311);
for MC = 1:MCLength
    P1 = plot(t_plot,std_psi_e_Con_MC(:,MC)*r2d*SM,'b');hold on; grid on;
    plot(t_plot,-std_psi_e_Con_MC(:,MC)*r2d*SM,'b')
end
for MC = 1:MCLength
    P2 =plot(t_plot,psi_Err_Con_MC(:,MC)*r2d,'r'); 
end
set(gca,'FontSize',axes_font_size)
legend([P1 P2],{[num2str(SM),'\sigma'],' CON \delta\Psi_3_2_1'},'FontSize',leg_font_size)
ylabel('$$\delta$$$$\psi$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
ylim(YLimitEAngCon)
EulErrAx(2) = subplot(312);
for MC = 1:MCLength
    plot(t_plot,std_the_e_Con_MC(:,MC)*r2d*SM,'b');hold on; grid on;
    plot(t_plot,-std_the_e_Con_MC(:,MC)*r2d*SM,'b')
end
for MC = 1:MCLength
    plot(t_plot,the_Err_Con_MC(:,MC)*r2d,'r'); 
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$\theta$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
ylim(YLimitEAngCon)
EulErrAx(3) = subplot(313);
for MC = 1:MCLength
    plot(t_plot,std_phi_e_Con_MC(:,MC)*r2d*SM,'b');hold on; grid on;
    plot(t_plot,-std_phi_e_Con_MC(:,MC)*r2d*SM,'b')
end
for MC = 1:MCLength
    plot(t_plot,phi_Err_Con_MC(:,MC)*r2d,'r'); 
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$\phi$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(EulErrAx,'xy')
ylim(YLimitEAngCon)
xlim(XLimit)
fig=gcf;
saveas(fig,[SaveLoc,'\CONEulAngErr',num2str(MCLength),'runs','.png'])

% TAC Euler Angle Errors %
fig=gcf;
figure(fig.Number+1)
EulErrAx(1) = subplot(311);
for MC = 1:MCLength
    P1 = plot(t_plot,std_psi_e_Tac_MC(:,MC)*r2d*SM,'b');hold on; grid on;
    plot(t_plot,-std_psi_e_Tac_MC(:,MC)*r2d*SM,'b')
end
for MC = 1:MCLength
    P2 =plot(t_plot,psi_Err_Tac_MC(:,MC)*r2d,'r'); 
end
set(gca,'FontSize',axes_font_size)
legend([P1 P2],{[num2str(SM),'\sigma'],' TAC \delta\Psi_3_2_1'},'FontSize',leg_font_size)
ylabel('$$\delta$$$$\psi$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
ylim(YLimitEAngTacNav)
EulErrAx(2) = subplot(312);
for MC = 1:MCLength
    plot(t_plot,std_the_e_Tac_MC(:,MC)*r2d*SM,'b');hold on; grid on;
    plot(t_plot,-std_the_e_Tac_MC(:,MC)*r2d*SM,'b')
end
for MC = 1:MCLength
    plot(t_plot,the_Err_Tac_MC(:,MC)*r2d,'r'); 
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$\theta$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
ylim(YLimitEAngTacNav)
EulErrAx(3) = subplot(313);
for MC = 1:MCLength
    plot(t_plot,std_phi_e_Tac_MC(:,MC)*r2d*SM,'b');hold on; grid on;
    plot(t_plot,-std_phi_e_Tac_MC(:,MC)*r2d*SM,'b')
end
for MC = 1:MCLength
    plot(t_plot,phi_Err_Tac_MC(:,MC)*r2d,'r'); 
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$\phi$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(EulErrAx,'xy')
ylim(YLimitEAngTacNav)
xlim(XLimit)
fig=gcf;
saveas(fig,[SaveLoc,'\TACEulAngErr',num2str(MCLength),'runs','.png'])

% NAV Euler Angle Errors %
fig=gcf;
figure(fig.Number+1)
EulErrAx(1) = subplot(311);
for MC = 1:MCLength
    P1 = plot(t_plot,std_psi_e_Nav_MC(:,MC)*r2d*SM,'b');hold on; grid on;
    plot(t_plot,-std_psi_e_Nav_MC(:,MC)*r2d*SM,'b')
end
for MC = 1:MCLength
    P2 =plot(t_plot,psi_Err_Nav_MC(:,MC)*r2d,'r'); 
end
set(gca,'FontSize',axes_font_size)
legend([P1 P2],{[num2str(SM),'\sigma'],' NAV \delta\Psi_3_2_1'},'FontSize',leg_font_size)
ylabel('$$\delta$$$$\psi$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
ylim(YLimitEAngTacNav)
EulErrAx(2) = subplot(312);
for MC = 1:MCLength
    plot(t_plot,std_the_e_Nav_MC(:,MC)*r2d*SM,'b');hold on; grid on;
    plot(t_plot,-std_the_e_Nav_MC(:,MC)*r2d*SM,'b')
end
for MC = 1:MCLength
    plot(t_plot,the_Err_Nav_MC(:,MC)*r2d,'r'); 
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$\theta$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
ylim(YLimitEAngTacNav)
EulErrAx(3) = subplot(313);
for MC = 1:MCLength
    plot(t_plot,std_phi_e_Nav_MC(:,MC)*r2d*SM,'b');hold on; grid on;
    plot(t_plot,-std_phi_e_Nav_MC(:,MC)*r2d*SM,'b')
end
for MC = 1:MCLength
    plot(t_plot,phi_Err_Nav_MC(:,MC)*r2d,'r'); 
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$\phi$$ (deg)','Interpreter','Latex','FontSize',yfont_size);
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(EulErrAx,'xy')
ylim(YLimitEAngTacNav)
xlim(XLimit)
fig=gcf;
saveas(fig,[SaveLoc,'\NAVEulAngErr',num2str(MCLength),'runs','.png'])

% CON Quaternion Errors %
fig=gcf;
figure(fig.Number+1)
QErrAx(1) = subplot(311);
for MC=1:MCLength
    P1 = plot(t_plot,q_s_s_hat1_Con_MC(:,MC),'r'); hold on; grid on;
    P2 = plot(t_plot,SM*std_qe1_Con_MC(:,MC),'b');
    plot(t_plot,-SM*std_qe1_Con_MC(:,MC),'b')
end
set(gca,'FontSize',axes_font_size)
legend([P1 P2],{[num2str(SM),'\sigma'],' CON \deltaq'},'FontSize',leg_font_size)
ylabel('$$\delta$$$$q_{S\hat{S},1}$$','Interpreter','Latex','FontSize',font_size);
QErrAx(2) = subplot(312);
for MC=1:MCLength
    P1 = plot(t_plot,q_s_s_hat2_Con_MC(:,MC),'r'); hold on; grid on;
    P2 = plot(t_plot,SM*std_qe2_Con_MC(:,MC),'b');
    plot(t_plot,-SM*std_qe2_Con_MC(:,MC),'b')
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$q_{S\hat{S},2}$$','Interpreter','Latex','FontSize',font_size);
QErrAx(3) = subplot(313);
for MC=1:MCLength
    P1 = plot(t_plot,q_s_s_hat3_Con_MC(:,MC),'r'); hold on; grid on;
    P2 = plot(t_plot,SM*std_qe3_Con_MC(:,MC),'b');
    plot(t_plot,-SM*std_qe3_Con_MC(:,MC),'b')
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$q_{S\hat{S},3}$$','Interpreter','Latex','FontSize',font_size);
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(QErrAx,'xy')
xlim(XLimit)
ylim([-0.4 0.4])
fig=gcf;
saveas(fig,[SaveLoc,'\CONQuatErr',num2str(MCLength),'runs','.png'])

% TAC Quaternion Errors %
fig=gcf;
figure(fig.Number+1)
QErrAx(1) = subplot(311);
for MC=1:MCLength
    P1 = plot(t_plot,q_s_s_hat1_Tac_MC(:,MC),'r'); hold on; grid on;
    P2 = plot(t_plot,SM*std_qe1_Tac_MC(:,MC),'b');
    plot(t_plot,-SM*std_qe1_Tac_MC(:,MC),'b')
end
set(gca,'FontSize',axes_font_size)
legend([P1 P2],{[num2str(SM),'\sigma'],' TAC \deltaq'},'FontSize',leg_font_size)
ylabel('$$\delta$$$$q_{S\hat{S},1}$$','Interpreter','Latex','FontSize',font_size);
QErrAx(2) = subplot(312);
for MC=1:MCLength
    P1 = plot(t_plot,q_s_s_hat2_Tac_MC(:,MC),'r'); hold on; grid on;
    P2 = plot(t_plot,SM*std_qe2_Tac_MC(:,MC),'b');
    plot(t_plot,-SM*std_qe2_Tac_MC(:,MC),'b')
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$q_{S\hat{S},2}$$','Interpreter','Latex','FontSize',font_size);
QErrAx(3) = subplot(313);
for MC=1:MCLength
    P1 = plot(t_plot,q_s_s_hat3_Tac_MC(:,MC),'r'); hold on; grid on;
    P2 = plot(t_plot,SM*std_qe3_Tac_MC(:,MC),'b');
    plot(t_plot,-SM*std_qe3_Tac_MC(:,MC),'b')
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$q_{S\hat{S},3}$$','Interpreter','Latex','FontSize',font_size);
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(QErrAx,'xy')
xlim(XLimit)
ylim([-0.1 0.1])
fig=gcf;
saveas(fig,[SaveLoc,'\TACQuatErr',num2str(MCLength),'runs','.png'])

% NAV Quaternion Errors %
fig=gcf;
figure(fig.Number+1)
QErrAx(1) = subplot(311);
for MC=1:MCLength
    P1 = plot(t_plot,q_s_s_hat1_Nav_MC(:,MC),'r'); hold on; grid on;
    P2 = plot(t_plot,SM*std_qe1_Nav_MC(:,MC),'b');
    plot(t_plot,-SM*std_qe1_Nav_MC(:,MC),'b')
end
set(gca,'FontSize',axes_font_size)
legend([P1 P2],{[num2str(SM),'\sigma'],' NAV \deltaq'},'FontSize',leg_font_size)
ylabel('$$\delta$$$$q_{S\hat{S},1}$$','Interpreter','Latex','FontSize',font_size);
QErrAx(2) = subplot(312);
for MC=1:MCLength
    P1 = plot(t_plot,q_s_s_hat2_Nav_MC(:,MC),'r'); hold on; grid on;
    P2 = plot(t_plot,SM*std_qe2_Nav_MC(:,MC),'b');
    plot(t_plot,-SM*std_qe2_Nav_MC(:,MC),'b')
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$q_{S\hat{S},2}$$','Interpreter','Latex','FontSize',font_size);
QErrAx(3) = subplot(313);
for MC=1:MCLength
    P1 = plot(t_plot,q_s_s_hat3_Nav_MC(:,MC),'r'); hold on; grid on;
    P2 = plot(t_plot,SM*std_qe3_Nav_MC(:,MC),'b');
    plot(t_plot,-SM*std_qe3_Nav_MC(:,MC),'b')
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$q_{S\hat{S},3}$$','Interpreter','Latex','FontSize',font_size);
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(QErrAx,'xy')
xlim(XLimit)
ylim([-0.1 0.1])
fig=gcf;
saveas(fig,[SaveLoc,'\NAVQuatErr',num2str(MCLength),'runs','.png'])


% Angular Error Between S3 and S3_hat, CON %
fig=gcf;
figure(fig.Number+1)
for MC=1:MCLength
    P2 = plot(t_plot,epsilon_Con_MC(:,MC),'r');hold on;
end
P1 = plot(t_plot,std_epsilon_Con_MC(:,end)*SM,'b','LineWidth',1);
set(gca,'FontSize',axes_font_size)
grid on; ylabel('$$\epsilon$$ (deg)','Interpreter','Latex','FontSize',font_size);
legend([P1 P2],{[num2str(SM),'\sigma'],'CON \epsilon'},'FontSize',leg_font_size)
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
xlim(XLimit)
ylim([0 20])
fig=gcf;
saveas(fig,[SaveLoc,'\CONEpsilon',num2str(MCLength),'runs','.png'])

% Angular Error Between S3 and S3_hat, TAC %
fig=gcf;
figure(fig.Number+1)
for MC=1:MCLength
    P2 = plot(t_plot,epsilon_Tac_MC(:,MC),'r');hold on;
end
P1 = plot(t_plot,std_epsilon_Tac_MC(:,end)*SM,'b','LineWidth',1);
set(gca,'FontSize',axes_font_size)
grid on; ylabel('$$\epsilon$$ (deg)','Interpreter','Latex','FontSize',font_size);
legend([P1 P2],{[num2str(SM),'\sigma'],'TAC \epsilon'},'FontSize',leg_font_size)
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
xlim(XLimit)
ylim([0 2])
fig=gcf;
saveas(fig,[SaveLoc,'\TACEpsilon',num2str(MCLength),'runs','.png'])

% Angular Error Between S3 and S3_hat, NAV %
fig=gcf;
figure(fig.Number+1)
for MC=1:MCLength
    P2 = plot(t_plot,epsilon_Nav_MC(:,MC),'r');hold on;
end
P1 = plot(t_plot,std_epsilon_Nav_MC(:,end)*SM,'b','LineWidth',1);
set(gca,'FontSize',axes_font_size)
grid on; ylabel('$$\epsilon$$ (deg)','Interpreter','Latex','FontSize',font_size);
legend([P1 P2],{[num2str(SM),'\sigma'],'NAV \epsilon'},'FontSize',leg_font_size)
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
xlim(XLimit)
ylim([0 2])
fig=gcf;
saveas(fig,[SaveLoc,'\NAVEpsilon',num2str(MCLength),'runs','.png'])



