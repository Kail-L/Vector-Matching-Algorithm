%% Constants file %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File to define relevant constants and parameters for the simulation.
%
% Kail Laughlin
% 7/20/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Unit Conversions %%
d2r=pi/180;     % degrees to radians conversion
r2d=1/d2r;      % radians to degrees conversion
g2mg = 1e3;     % gauss to milligauss conversion
mg2g = 1/g2mg;  % milligauss to gauss conversion
T2mg = 1e7;     % tesla to milligauss conversion
mg2T = 1/T2mg;  % milligauss to tesla conversion
T2g = 1e4;      % tesla to gauss conversion
g2T = 1/T2g;    % gauss to tesla conversion
s2hr = 1/3600;  % seconds to hours conversion
hr2s =  1/s2hr; % hours to seconds conversion

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Orbital Parameters %%
const.mu1 = 3.986e14;       % Earth gravitational constant (m^3/s^2)
const.J2 = 1.08262645e-3;   % J2 perturbation coefficient (set to 0 if no J2 perturbation)
const.Re = 6371.2e3;        % radius of the Earth (m)
const.a =const.Re+408000;   % semimajor axis of orbit (m)
const.e = 0.001;            % eccentricity of orbit
const.inc =50*pi/180;       % inclination of orbit (rad)
const.omega = 0;            % argument of perigeee (rad)
const.Omega = 0;            % right ascension of the ascending node (rad)
const.t0 = 0;               % time of perigee passage (s)
const.T = round(2*pi/...
          sqrt(const.mu1)*const.a^(3/2));  % Orbital Period (s)
% NOTE: The orbital period calculated is for an elliptical orbit. %

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Spacecraft Parameters %%
const.m_bscal = 0;%1e-4;       % residual magnetic dipole moment of spacecraft (A*m^2)
const.m_b = const.m_bscal*[1;1;0];    
const.UseSGPB = 0;             % boolean to use SGPB data (Need to be UMN student to obtain)
const.UseHyCUBE = 0;           % boolean to use HyCUBE reentry data (note no mag data included, only inertial AD)
if const.UseSGPB ==1 || const.UseHyCUBE ==1
    const.UseExternal = 1;     % boolean to use external sim data seen in ExternalData.m
else
    const.UseExternal = 0;     % boolean to use external sim data seen in ExternalData.m
end
const.SatLength = 10;          % length of cubesat, cm
const.SatWidth = 10;           % Width of cubesat, cm
const.SatDepth = 30;           % Depth of cubesat, cm

if const.UseSGPB == 0          % SOCRATES Heritage MOI info
    const.m_s = 3.330;         % mass of satellite (kg)
    Ixx = 0.0508;   % (kg m^2)
    Iyy = 0.0508;   % (kg m^2)
    Izz = 0.0252;   % (kg m^2)
    Ixy = 0;        % (kg m^2)
    Iyz = 0;        % (kg m^2)
    Ixz = 0;        % (kg m^2)
else                           % SGPB MOI info
    const.m_s = 3154.82;       % mass of satellite (kg)
    Ixx = 5138;     % (kg m^2)
    Iyy = 5065;     % (kg m^2)
    Izz = 3450;     % (kg m^2)
    Ixy = -4.23;    % (kg m^2)
    Iyz = 0;        % (kg m^2)
    Ixz = 0;        % (kg m^2)
end     
const.I_b = [Ixx Ixy Ixz
             Ixy Iyy Iyz
             Ixz Iyz Izz]; % (kg m^2)
         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sim Initialization Parameters %%
const.NumOrbits = 5;         % Number of orbits to simulate
const.MCLength = 1;          % Number of Monte Carlo runs to perform
const.UseCustomStartDate = 0;% Boolean to decide if custom sim start date given below used
const.SimStartDate = [2019 5 22 0 0 1]; % Simulation Start Date [year month day hour minute sec]

% Daily Reboot %
% Note: Daily reboot is meant to simulate the planned satellite reset that
% occurs every 24 hours for EXACT & IMPRESS. 
const.PerformDailyReboot = 0;% Boolean to decide if daily reboot performed
const.DailyRebootTime = 6;   % When, in hours, the daily reboot occurs from initialization
const.RebootLength = 2*60;   % Length of reboot in seconds

% Magnetic Field Model %
% The magnetic field model used is dependent on the function used. A less
% accurate but faster model is EarthMagField.m (adapted from Richard Forbes
% and Ryan Caverly), which takes position in ECI frame as an input. The
% world magnetic model (WMM) is an international model developed by the UK
% and is a professional MATLAB function. The WMM however takes much longer
% times in terms of computation, and takes position (in LLA coords) as an
% input. The international geomagnetic reference model (IGRF) is another
% quality model, but requires up to date info that MATLAB may not have.
%
% 0 - EarthMagField.m is used
% 1 - World Magnetic Model (WMM) is used (NOTE: MATLAB 2020a is required)
% 2 - International Geomagnetic Reference Field (IGRF) is used (NOTE: MATLAB 2020a is required)

const.MagFieldFlag = 1;      

% Initial Angular Rates (rad/s) %

% const.omega_S_ECI_0 = [0;0;0];               % Zero angular velocity
% const.omega_S_ECI_0 = [0;0;0.001];           % Low spin rate
% const.omega_S_ECI_0 = [0;0;0.01];            % Spin stabilized rate
const.omega_S_ECI_0 = [0.02;-0.05;0.03];     % Large initial rates (Realistic)
% const.omega_S_ECI_0 = [0.3;0.3;0.3];         % Very high intial rates (limit?)

% Initial Truth Attitude %
const.Eul0=[0; 80; 0]*pi/180;         % Initial Truth 321 Euler Angles, [yaw pitch roll]^T in rad

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EKF Initialization Parameters %%
% const.IQ = 0 -> Uses true quaternion as initial quaternion estimate
%
% const.IQ = 1 -> Uses a "Maximum Angular Error" to create an error DCM
%                 that is used to rotate the true attitude to some error
%                 attitude. Error is bounded by the magnitude of const.ME.
%
% const.IQ = 2 -> Uses a completely random quaternion as initial estimate.
%
% const.IQ = 3 -> Uses a set of Error euler angles to initialize quaternion

const.IQ = 2;           % Flag for quaternion Initialization
const.ME = 90;          % Maximum angular error, deg, IQ = 1
const.ErrEul = [0; 160; 0]*pi/180;% Euler angle error

const.UseMeas = 1;      % Boolean to perform TAM measurement update
const.EstimateSF = 1;   % Boolean to estimate scale factors
const.constrainEst = 1; % Boolean to constrain scale factor estimates using QuadProg
const.RunSOTest = 1;    % Boolean to run Observability tests

% Initial Covariances %
if const.EstimateSF ==1
    const.P0 = [0.3^2*eye(3) zeros(3) zeros(3);
                zeros(3) (0.2*d2r)^2*eye(3) zeros(3);
                zeros(3) zeros(3) (0.02)^2*eye(3)];
            
else
    const.P0 = [0.3^2*eye(3) zeros(3);
                zeros(3) (0.2*d2r)^2*eye(3)];
end
% Tuning Parameters on G Matrix %
if const.UseSGPB ==1
    const.alpha = 1;    % Must be greater than 0. Set equal to 1 for no tuning. 1 works well.
    const.beta = 1;   % Must be greater than 0. Set equal to 1 for no tuning. 6-10 works well.
    const.gamma = 1;    % Must be greater than 0. Set equal to 1 for no tuning. 1 works well.
else
    const.alpha = 1;    % Must be greater than 0. Set equal to 1 for no tuning. 1 works well.
    const.beta = 1;     % Must be greater than 0. Set equal to 1 for no tuning. 6-10 works well.
    const.gamma = 1;    % Must be greater than 0. set equal to 1 for no tuning. 1 works well.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Sensor Parameters %%

% Sample Rates %
const.GyroSamp = 10;            % Sample rate of gyro, Hz
const.MagSamp = 0.1;            % Sample rate of magnetometer, Hz
const.SunSamp = 0.1;            % Sample rate of sun sensor, Hz
const.SampRatio = const.GyroSamp/const.MagSamp; % Time update to measurement update ratio

% Gyro Params %
% Note: const.G_SF and const.G_biasRepeat are bounds used to initialize
% random values of the scale factors and biases for each run of the VMA. 
const.G_AddNullShift = 1;       % Boolean to add nullshift 
const.G_AddWhiteNoise = 1;      % Boolean to add white noise
const.G_AddGaussMarkovNoise = 1;% Boolean to add Gauss Markov noise
const.G_AddScaleFactor = 1;     % Boolean to add Scale Factor
const.G_UseSetVals = 0;         % Boolean to emulate bias and scale factor with set values
const.G_UseARW = 1;             % Boolean to use ARW to create white noise on system and to use in Qw(t)
const.G_ScaleFactorLB = -0.05;  % Lower bound scale factor estimate for quadprog
const.G_ScaleFactorUB = 0.05;   % Upper bound on scale factor estimate for quadprog
const.G_tau = 300;              % Gyro gauss markov time constant, s 
const.G_ARW = 0.26;             % Gyro angle random walk deg/sqrt(hr)
const.G_inRunBiasSTD = 5.1;     % Gyro in-run bias stability, deg/hr
const.G_SF = 0.05;              % Gyro scale factor limit, percent
const.G_biasRepeat = 0.4;       % Gyro bias repeatability limit, deg/s
const.G_measNoise = 0.135;      % Gyro output noise, deg/s rms
const.G_SFNoise = 200;          % Gyro SF std, ppm
const.CSg = [1 0 0;             % Misalignment DCM between gyro and body 
             0 1 0;
             0 0 1];
if const.G_UseSetVals == 1
    const.G_SFVal = [0.001 0.0015 -0.002]; % Constant scale factors (percentages) (only if const.G_UseSetVals = 1)
    const.G_b0Val = [0.2 0.2 0.2];    % Constant b0 value (deg/s) (only if const.G_UseSetVals = 1)
end


% Mag Params %
const.M_AddHardIron = 0;        % Boolean to add hard iron biases
const.M_AddSoftIron = 0;        % Boolean to add soft iron biases
const.M_AddScaleFactor = 0;     % Boolean to add scale factor
const.M_AddWhiteNoise = 1;      % Boolean to add white noise
const.M_UseSetVals = 0;         % Boolean to emulate bias and scale factor with set values
const.M_ScaleFactorLB = 0.95;   % Lower bound scale factor estimates
const.M_ScaleFactorUB = 1.05;   % Upper bound on scale factor estimates
const.M_SF = 0.017;             % Gyro scale factor standard deviation, percent
const.M_HardIron = 20;         % Hard iron bias standard deviation, mgauss
const.M_alpha11 = 0.05;         % Soft iron induction parameter, 1,1 direction
const.M_alpha12 = 0.10;         % Soft iron induction parameter, 1,2 direction
const.M_alpha13 = 0.01;         % Soft iron induction parameter, 1,3 direction
const.M_alpha21 = 0.10;         % Soft iron induction parameter, 2,1 direction
const.M_alpha22 = 0.03;         % Soft iron induction parameter, 2,2 direction
const.M_alpha23 = 0.20;         % Soft iron induction parameter, 2,3 direction
const.M_alpha31 = 0.01;         % Soft iron induction parameter, 3,1 direction
const.M_alpha32 = 0.20;         % Soft iron induction parameter, 3,2 direction
const.M_alpha33 = 0.05;          % Soft iron induction parameter, 3,3 direction

if const.M_UseSetVals == 1
    const.M_SFVal = [0.1 0.05 -0.02]; % Constant scale factors (percentages) (only if const.G_UseSetVals = 1)
    const.M_HIVal = [0.2 0.2 0.2];    % Constant hard iron value (mGauss) (only if const.G_UseSetVals = 1)
end

if const.UseSGPB ==1
    const.M_measNoise = 1.2;   % Mag sensor variance for cov. matrix, mgauss, SGPB
else
    const.M_measNoise = 1.1;    % Mag sensor variance for cov. matrix, mgauss  
end

% Sun Params %
const.UseSun = 0;               % Boolean to use Sun Sensor
const.S_measNoise = 6;          % Error in sun vector measurement, deg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plotting & Visualizer Parameters %%
direc = pwd;
if const.UseSun == 1
    const.VecNum=2;
else
    const.VecNum=1;
end
const.plotMean=0; %Boolean to plot Monte Carlo Mean
const.plotTime=1; %1 for hour, 0 for seconds
const.KnowledgeError = 10;      % Used in Angular Error Plot
const.SigmaMultiplier = 3;      % Used in any plot with covariance bound
const.SaveLocation = ['C:\Users\'];
const.CreateVid = 0;            % Boolean to create video of cubesat
const.VidTime = 10;             % Length of video, seconds
