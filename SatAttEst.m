%% Satellite Att. Estimation Code %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SatAttEst.m provides "Truth" attitude data on a satellite with given 
% parameters then corrupts it based on sensor information. An EKF using 
% single vector(magnetometer) measurement OR dual vector (magnetometer and
% sun vector) measurement is used to estimate the attitude of a satellite.
%
% To change or modify parameters of simulation, update values within
% constants_struct.m. This struct inlcudes orbit, spacecraft, and sensor
% parameters.
%
% Kail Laughlin, 7/7/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Function Flow %%
% Based upon the relevant parameters within constants_struct.m, This script
% pushes data through ODE45 to integrate satellite attitude and
% orbital info for a given time frame. This integrated data is considered
% the "Truth" data that is then corrupted and used within the EKF
% algorithm. A flow of the highlevel functionality of the sim is given
% below:
%
% constants_struct (line 65)
%   - Call the relevant constants to initialize the sim
%
% Initialization (lines 67 - 300)
%   -Initialize info for simulation and for display to command window
%
% Integration (lines 303 - 314)
%   - Use ODE45 to integrate dynamic equations set up in ODEs.m (if
%     external data is not used)
%
% Post Processing (lines 319 - 329)
%   - Take date from ODE45 and package into useable format
%
% Extended Kalman Filter (lines 365 - 373)
%   - Depending on choice of EKF algorithm, a 1 vector or 2 vector solution
%     will be used
%   - Corrupts data via emulateGyro.m and emulateMag.m near line 177
%   - Uses VMA_EKF.m or VMA_EKF_MC.m depending on whether Monte Carlo runs
%     are being performed
%   - If the observability test is run, it is called on line 650 of
%     VMA_EKF.m and line 660 of VMA_EKF_MC.m
%
% Visualizer (line 379)
%   - Will create a *.avi movie of the cubesat rotating in inertial frame
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear 
close all
format short

% Font size, line size, and line width for plotting %
font_size = 10; 
line_size = 10;
line_width = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add subfolders to filepath %%
p = pwd;
folder = genpath(p);
addpath(folder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call Constants Structure %%
constants_struct    % All constants in one file. 

%% Use External Data? %%
if const.UseExternal == 1   % If using external data, don't use ODE45 output
    disp('Pre Processing External Data...')
    Total_Time = 0;
    tic
    ExternalData                % Run pre processing script
    ExternalDataTime = toc;
    Total_Time= Total_Time+ExternalDataTime;
    disp(['Runtime: ',num2str(Total_Time),' s'])
    q_S_ECI_0 = x_out(1,7:10)';    
    C_S_ECI_0 = eul2dcm(q_S_ECI_0);
    eul = dcm2eul(C_S_ECI_0);
    om_S_ECI_0 = x_out(1,11:13);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Sim Info to Command Window %   
    disp('------------------Sim Parameters----------------')
    disp('Using External Data')
    if const.MCLength == 1
    disp('Run Monte Carlo?:           No')
    else
    disp('Run Monte Carlo?:           Yes')
    disp(['Number of Monte Carlo Runs: ',num2str(const.MCLength)])
    end
    fprintf('\n')
    disp('---------------Sensor Parameters----------------')
    if const.UseMeas == 1
    disp(['Mag Sample Rate:            ',num2str(const.MagSamp),' Hz'])
    disp(['Mag Output Noise:           ',num2str(const.M_measNoise),' mGauss'])
    end
    disp(['Gyro Sample Rate:           ',num2str(const.GyroSamp),' Hz'])
    if const.G_AddGaussMarkovNoise ==1
    disp(['Gyro In-Run Bias:           ',num2str(const.G_inRunBiasSTD),' deg/hr'])
    disp(['Gyro Correlation Time:      ',num2str(const.G_tau),' s'])
    end
    if const.G_AddNullShift ==1
    disp(['Gyro Constant Bias Limit:   ',num2str(const.G_biasRepeat),' deg/s'])
    end
    if const.G_AddScaleFactor ==1
    disp(['Gyro Scale Factor Limit:    ',num2str(const.G_SF),])
    end
    if const.G_AddWhiteNoise ==1
    disp(['Gyro ARW:                   ',num2str(const.G_ARW),' deg/sqrt(hr)'])
    disp(['Gyro Output Noise:          ',num2str(const.G_measNoise),' deg/s'])
    end
    if const.UseSun ==1
    disp(['Sun Sens Sample Rate:       ',num2str(const.SunSamp),' Hz'])
    disp(['Sun Sensor Meas Error:      ',num2str(const.S_measNoise),' deg'])
    end
    fprintf('\n')
    disp('---------------Attitude Parameters--------------')
    disp(['Initial Yaw:                ',num2str(eul(1)*180/pi),' deg'])
    disp(['Initial Pitch:              ',num2str(eul(2)*180/pi),' deg'])
    disp(['Initial Roll:               ',num2str(eul(3)*180/pi),' deg'])
    disp(['Initial p (about 1 axis):   ',num2str(om_S_ECI_0(1)*180/pi),' deg/s'])
    disp(['Initial q (about 2 axis):   ',num2str(om_S_ECI_0(2)*180/pi),' deg/s'])
    disp(['Initial r (about 3 axis):   ',num2str(om_S_ECI_0(3)*180/pi),' deg/s'])
    fprintf('\n')
    disp('-----------------EKF Parameters-----------------')
    if const.UseSun ==1
    disp('EKF Sensors:                IMU, TAM, Sun Sensor') 
    else
    disp('EKF Sensors:                IMU, TAM')
    end
    if const.G_UseARW ==0
    disp('Driving noise in Qw(t):     Output noise, deg/s')
    else
    disp('Driving noise in Qw(t):     ARW, deg/sqrt(s)')
    end
    if const.RunSOTest ==1
    disp('Run Observability Test?:    Yes')
    else
    disp('Run Observability Test?:    No')
    end
    fprintf('\n')
    fprintf('\n')
else
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Orbital Conditions %%
% Determine starting position (m) and velocity (m/s) %
    [r_a0,v_a0] = findRandV(const.mu1,const.e,const.inc,...
                            const.a,const.omega,const.Omega,const.t0,0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial Attitude Conditions %%
% Can initialize true att. with random quaternion or fixed euler angles %
% InitAtt = 0 -> Initialize with set Euler Angles (321 rotation sequence)
% InitAtt = 1 -> Initialize with random quaternion
    InitAtt = 0;              

    if InitAtt ==0
        eul = const.Eul0;     % [yaw pitch roll]^T in rad
        C_S_ECI_0 = eul2dcm(eul);
        q_S_ECI_0 = eul2quat(eul);
    else                    % Random quat
        q0 = rand(1,1);
        q1 = 2*rand(1,1) - 1;
        q2 = 2*rand(1,1) - 1;
        q3 = 2*rand(1,1) - 1;
        q_S_ECI_0 = [q0 q1 q2 q3]';
        q_S_ECI_0 = q_S_ECI_0./norm(q_S_ECI_0);

        C_S_ECI_0=quat2dcm(q_S_ECI_0');
        eul = dcm2eul(C_S_ECI_0);
    end
    
    q_0 = q_S_ECI_0(1);                % Quaternion Scalar
    q_bar = q_S_ECI_0(2:4);            % Quaternion Vector

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Initial Rotational Conditions %%
    om_S_ECI_0 = const.omega_S_ECI_0;      % Initial angular rates of sat
    h_S_ECI_0 = const.I_b*om_S_ECI_0;         % Initial angular momentum of sat

    % Initial condition vector for ODE45 % 
    IC = [r_a0; v_a0; q_0; q_bar; om_S_ECI_0; h_S_ECI_0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation Time & Timing Vectors %%
    t0 = 0; % s
    t_max = const.T*const.NumOrbits; % s 

    % Sensor timing vectors %
    t_span = (t0:1/const.GyroSamp:t_max);
    t_mag = (t0:1/const.MagSamp:t_max)';
    t_sun = (t0:1/const.SunSamp:t_max)';

    % Get UTC time of sim %
    c = clock;    % c = [year month day hour minute seconds] all integers

    if const.UseCustomStartDate == 1
        SimStart = datetime(const.SimStartDate);
    else
        % Initialize timing at 12:01 AM of current day %
        SimStart = datetime([c(1) c(2) c(3) 0 0 1]);
    end

    t_utc = datevec(SimStart+seconds(t_span)); 
    t_utc_mag = datevec(SimStart+seconds(t_mag));
    t_utc_sun = datevec(SimStart+seconds(t_sun)); 

    % Initial Geodetic (latitude, longitude, altitude) Coordinates %
    LLA_0 = eci2lla(r_a0',t_utc(1,:));
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Display Sim Info to Command Window %
    disp('------------------Sim Parameters----------------')
    disp(['Simulation Initial Time:    ',datestr(SimStart)])
    disp(['Simulation Timespan:        ',num2str(t_max),' s = ',num2str(t_max/3600),' hrs'])
    disp(['Number of Time Steps:       ',num2str(length(t_span))])
    disp(['Length of Gyro Time Vector: ',num2str(length(t_span))])
    if const.UseMeas == 1
    disp(['Length of Mag Time Vector:  ',num2str(length(t_mag))])
    end
    if const.UseSun == 1
    disp(['Length of Sun Time Vec:     ',num2str(length(t_sun))])
    end
    if const.MCLength == 1
    disp('Run Monte Carlo?:           No')
    else
    disp('Run Monte Carlo?:           Yes')
    disp(['Number of Monte Carlo Runs: ',num2str(const.MCLength)])
    end
    fprintf('\n')
    disp('---------------Spacecraft Parameters------------')
    disp(['Residual Magnetic Dipole:   ',num2str(const.m_bscal),' A*m^2'])
    disp(['Mass:                       ',num2str(const.m_s),' kg'])
    disp('Moment of Inertia Matrix (kg*m^2):')
    I_b=const.I_b
    disp('---------------Orbital Parameters---------------')
    disp(['Initial Latitude:           ',num2str(LLA_0(1)),' deg'])
    disp(['Initial Longitude:          ',num2str(LLA_0(2)),' deg'])
    disp(['Initial Altitude:           ',num2str((r_a0(1,1)-const.Re)/1000),' km'])
    disp(['Initial Velocity:           [',num2str(v_a0(1)/1000),' ',num2str(v_a0(2)/1000),' ',num2str(v_a0(3)/1000),'] km/s'])
    disp(['Inclination:                ',num2str(const.inc*180/pi),' deg'])
    disp(['Eccentricity:               ',num2str(const.e)])
    disp(['Orbital Period:             ',num2str(const.T/60),' min'])
    disp(['J2 Perturbation:            ',num2str(const.J2)])
    disp(['Argument of Perigee:        ',num2str(const.omega*180/pi),' deg'])
    disp(['Right Ascension:            ',num2str(const.Omega*180/pi),' deg'])
    disp(['Time of Perigee Passage:    ',num2str(const.t0),' s'])
    fprintf('\n')
    disp('---------------Attitude Parameters--------------')
    disp(['Initial Yaw:                ',num2str(eul(1)*180/pi),' deg'])
    disp(['Initial Pitch:              ',num2str(eul(2)*180/pi),' deg'])
    disp(['Initial Roll:               ',num2str(eul(3)*180/pi),' deg'])
    disp(['Initial p (about 1 axis):   ',num2str(om_S_ECI_0(1)*180/pi),' deg/s'])
    disp(['Initial q (about 2 axis):   ',num2str(om_S_ECI_0(2)*180/pi),' deg/s'])
    disp(['Initial r (about 3 axis):   ',num2str(om_S_ECI_0(3)*180/pi),' deg/s'])
    fprintf('\n')
    disp('---------------Sensor Parameters----------------')
    if const.UseMeas == 1
    disp(['Mag Sample Rate:            ',num2str(const.MagSamp),' Hz'])
    disp(['Mag Output Noise:           ',num2str(const.M_measNoise),' mGauss'])
    end
    disp(['Gyro Sample Rate:           ',num2str(const.GyroSamp),' Hz'])
    if const.G_AddGaussMarkovNoise ==1
    disp(['Gyro In-Run Bias:           ',num2str(const.G_inRunBiasSTD),' deg/hr'])
    disp(['Gyro Correlation Time:      ',num2str(const.G_tau),' s'])
    end
    if const.G_AddNullShift ==1
    disp(['Gyro Constant Bias Limit:   ',num2str(const.G_biasRepeat),' deg/s'])
    end
    if const.G_AddScaleFactor ==1
    disp(['Gyro Scale Factor Limit:    ',num2str(const.G_SF),])
    end
    if const.G_AddWhiteNoise ==1
    disp(['Gyro ARW:                   ',num2str(const.G_ARW),' deg/sqrt(hr)'])
    disp(['Gyro Output Noise:          ',num2str(const.G_measNoise),' deg/s'])
    end
    if const.UseSun ==1
    disp(['Sun Sens Sample Rate:       ',num2str(const.SunSamp),' Hz'])
    disp(['Sun Sensor Meas Error:      ',num2str(const.S_measNoise),' deg'])
    end
    fprintf('\n')
    disp('-----------------EKF Parameters-----------------')
    if const.UseSun ==1
    disp('EKF Sensors:                IMU, TAM, Sun Sensor') 
    else
    disp('EKF Sensors:                IMU, TAM')
    end
    if const.G_UseARW ==0
    disp('Driving noise in Qw(t):     Output noise, deg/s')
    else
    disp('Driving noise in Qw(t):     ARW, deg/sqrt(s)')
    end
    if const.RunSOTest ==1
    disp('Run Observability Test?:    Yes')
    else
    disp('Run Observability Test?:    No')
    end
    fprintf('\n')
    fprintf('\n')
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulation %%
    Total_Time = 0;
    options = odeset('AbsTol',1e-12,'RelTol',1e-12); % This changes the integration tolerence.
    disp('Creating Truth Data...')
    tic
    [t,x_out] = ode45(@ODEs,t_span,IC,options,const);% Integrate "truth" data
    Integration_Time = toc;
    Total_Time = Total_Time+Integration_Time;
    if Integration_Time<60
        disp(['Runtime: ',num2str(Integration_Time),' s'])
    else
        disp(['Runtime: ',num2str(Integration_Time/60),' min'])
    end
    fprintf('\n')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Post Processing %%
disp('Post Processing Data...')
tic
post_processing 
PostProcessing_Time = toc;
Total_Time = Total_Time+PostProcessing_Time;
if PostProcessing_Time<60
    disp(['Runtime: ',num2str(PostProcessing_Time),' s'])
else
    disp(['Runtime: ',num2str(PostProcessing_Time/60),' min'])
end
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot ODE45 Results
% PS = 0 -> No ODE45 results plotted
% PS = 1 -> ODE45 results plotted
%
% NOTE: Graphs that will be plotted include:
%      - 3D Spacecraft Orbit
%      - True position and velocity in inertial ECI frame
%      - True Attitude Quaternion
%      - True Euler angles
%      - True Quaternion
%      - Total Energy check
%      - True Body Rotational Rates, p,q, & r.

PS = 0;             

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot data %%
if PS==1 && const.UseExternal ==0
    disp('Plotting Truth Data...')
    tic
    plot_script
    Plotting_Time = toc;
    Total_Time = Total_Time+Plotting_Time;
    if Plotting_Time<60
        disp(['Runtime: ',num2str(Plotting_Time),' s'])
    else
        disp(['Runtime: ',num2str(Plotting_Time/60),' min'])
    end
    fprintf('\n')
end

% MagCalTestScript
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EKF %%
disp('Initializing EKF...')
tic
if const.MCLength==1
    VMA_EKF
else
    VMA_EKF_MC
end
Estimation_Time = toc;
Total_Time = Total_Time+Estimation_Time;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Visualizer %%
% NOTE: This takes a while to run %
if const.CreateVid ==1
    l=const.SatLength; w=const.SatWidth; d=const.SatDepth;
    t_anim = const.VidTime;
    disp('Running Visualizer...')
    fprintf('\n')
    tic
    AttitudeAnim('CubeSatAttitude',q_ba,qhat_ba,l,w,d,t,t_anim)
    Vis_Time=toc;
    Total_Time = Total_Time+Vis_Time;
end
disp(['Total Sim Runtime: ',num2str(Total_Time/60),' mins'])
