function [Gyro_Corrupt,G_SF,B_true] = emulateGyro(t,Gyro_Good,const)
%EMULATEGYRO Generates a corrupted gyroscope vector based on gyro
%characteristics defined in constants_struct.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GYRO_CORRUPT = emulateGyro(t,Gyro_Good,const) takes in a "clean" or true
% gyroscope measurement and corrupts it with potentially 4 sources of bias:
%
% 1.) Constant null shift
% 2.) Wide band "white" noise
% 3.) Gauss-Markov correlated noise
% 4.) Scale factor errors
%
% The choice of biases used is dependent on a flag within the const
% structure input into the function. Note that all computation occurs in
% radians.
%
% SOURCES:
% Built off of code from Demoz Gebre-Egziabher, but largely created from
% scratch. For some detailed info on the model used, see Laughlin
% "Single-Vector Aiding of an IMU for CubeSat Attitude Determination"
% Chapter 3.
%
%           C_Sg*W_corrupt = (1+SF)*(W_good) + C_Sg*(b_g(t) + n_g(t))
%
% where b_g is the time varying gyro bias given by
%
%                           b_g(t) = b_0 + b_1(t)
%
% where b_0 is the constant bias null shift and b_1 is the time varying
% gauss markov process noise. The variable n_g is the additive wideband
% noise. The (1+SF) value is a 3x3 matrix of scale factor errors.
%
% INPUT PARAMETERS:
% t = nx1 time vector used for gyro data, seconds
% Gyro_Good = nx3 array of true or "clean" gyro measurements, in rad/s
% const = constants structure with relevant gyro parameters
%
% OUTPUT PARAMETERS: 
% Gyro_Corrupt = nx3 array of biased or "noisy" gyro measurements, in rad/s
%
% VARIABLES:
% NullShift = nx3 array of nullshift constants that are added to clean gyro
% WhiteNoise = nx3 array of random gaussian variables multiplied by Output
%              noise standard deviation, G_measNoise
% GyroCorrelated = nx3 array of correlated noise variables. Each instance
%                  is correlated with the previous instances
% GyroSF = Scale Factored Gyro outputs
%
% Kail Laughlin
% Updated 6/28/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define some constants %%
d2r = pi/180;                     % Degrees to radians 
r2d = 1/d2r;                      % Radians to degrees
s2hr = 1/3600;                    % Seconds to hours
hr2s = 3600;                      % Hours to seconds
n = length(Gyro_Good);            % Length of gyro vector
dt = mean(diff(t));               % Time between gyro measurements

% Pull constants from constants_struct.m %
G_ARW = const.G_ARW;              % Gyro angular random walk (ARW) in deg/sqrt(hr)
G_IRb = const.G_inRunBiasSTD;     % Gyro In-Run bias repeatability, deg/hr
G_null = const.G_biasRepeat;      % Gyro bias repeatability, deg/s
G_tau = const.G_tau;              % Gyro correlation time for gauss-markov process
G_OutNoise = const.G_measNoise;   % Gyro output/measurement noise, deg/s rms
G_ScaleFactor = const.G_SF;       % Scale factor limit
G_SFnoise = const.G_SFNoise;      % Scale Factor 1 sigma std 
AddNull = const.G_AddNullShift;     % Boolean to add nullshift 
AddWN = const.G_AddWhiteNoise;      % Boolean to add white noise
AddGM = const.G_AddGaussMarkovNoise;% Boolean to add Gauss Markov noise
AddSF = const.G_AddScaleFactor;     % Boolean to add Scale Factor
UseSetVals = const.G_UseSetVals;    % Boolean to use set scale factors and b0 from constants structure
C_S_g = const.CSg;                      % DCM between gyro frame and spacecraft frame

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Constant Null Shift Errors %%
if UseSetVals==1
    NullShift = AddNull*(const.G_b0Val'*d2r);
    B0 = NullShift*r2d
    NullShift = (diag(NullShift)*ones(3,n))';
else
    NullShift = AddNull*(G_null*(2*rand(1,3)-1)*(d2r));
    B0 = NullShift*r2d
    NullShift = (diag(NullShift)*ones(3,n))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Wideband "White" Noise %%
if const.G_UseARW == 0
    WhiteNoise = AddWN*(G_OutNoise*(d2r)*randn(n,3));
else
%     WhiteNoise = AddWN*((G_ARW/60*d2r)*1/sqrt(330)*randn(n,3));
    WhiteNoise = AddWN*(G_ARW/60*d2r)*randn(n,3);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Gauss-Markov Correlated Noise %%
% Set up a state-space representation of noise %
a_gyro = -1/G_tau;
b_gyro = 1;
c_gyro = 1;
d_gyro = 0;

Q = 2*(G_IRb*(d2r)*(s2hr))^2/G_tau;  % Driving Noise White Power Spectral Density
Qd = disrw(a_gyro,b_gyro,dt,Q);      % Discretize spectral density

% State space of gyro %
SS_gyro = ss(a_gyro,b_gyro,c_gyro,d_gyro);

% Discrete time state space %
SS_gyro_dis = c2d(SS_gyro,dt);

% Discrete time A, B, C, & D matrices %
[ad,~,~,~] = ssdata(SS_gyro_dis);

sigmaU = sqrt(Qd);
ug = sigmaU*randn(length(t),3);
GyroCorrelated = zeros(length(ug),3);
adg = ad*eye(3);

%% Create Scale Factors %%
if UseSetVals == 1
    G_SF = const.G_SFVal';
    if AddSF==1
        G_SF
    end
else
    G_SF = G_ScaleFactor*(2*rand(1,3)-1);
    if AddSF==1
        G_SF
    end
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correlate Noise %%
for k=2:n
    GyroCorrelated(k,:) = AddGM*(adg*GyroCorrelated(k-1,:)' + ug(k-1,:)')';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Combine All Gyro Errors %%
Gyro_Corrupt = (diag(1+AddSF*G_SF)*Gyro_Good')' + (C_S_g*(GyroCorrelated + WhiteNoise + NullShift)')';
B_true = GyroCorrelated + WhiteNoise + NullShift;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end