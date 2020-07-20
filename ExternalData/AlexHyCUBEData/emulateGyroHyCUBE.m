function [Gyro_Corrupt,G_SF] = emulateGyroHyCUBE(t,Gyro_Good,IRb,Tau,OutNoise)
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
% scratch. For some detailed info on the model used, see Bageshwar et. al
% 2008, "Inertialy-Aided Vector Matching Algorithm for Attitude
% Determination of Spin Stabilized Satellites" pg. 5. We have added the
% ability to apply scale factors, so the bias model looks as follows
%
%           W_corrupt = (1+SF)*(W_good) + b_g(t) + n_g(t)
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
% Gyro_Good = nx3 array of true or "clean" gyro measurements, in rad
% const = constants structure with relevant gyro parameters
%
% OUTPUT PARAMETERS: 
% Gyro_Corrupt = nx3 array of biased or "noisy" gyro measurements, in rad
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
% Updated 4/3/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define some constants %%
d2r = pi/180;                     % Degrees to radians 
r2d = 1/d2r;                      % Radians to degrees
s2hr = 1/3600;                    % Seconds to hours
hr2s = 3600;                      % Hours to seconds
n = length(Gyro_Good);            % Length of gyro vector
dt = mean(diff(t));               % Time between gyro measurements

% Pull constants from constants_struct.m %

G_IRb = IRb;                      % Gyro In-Run bias repeatability, rad/s
G_null = 0;                       % Gyro bias repeatability, rad/s
G_tau = Tau;                        % Gyro correlation time for gauss-markov process
G_OutNoise = OutNoise;              % Gyro output/measurement noise, rad/s rms
G_SF = [0;0;0];                     % 3x1 array of gyro scale factors
AddNull = 0;                        % Boolean to add nullshift 
AddWN = 1;                          % Boolean to add white noise
AddGM = 1;                          % Boolean to add Gauss Markov noise
AddSF = 0;                          % Boolean to add Scale Factor

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Constant Null Shift Errors %%
% Randomize if Monte Carlo being run, use set value for single run.
NullShift = AddNull*((G_null*randn(1,3)));
NullShift = (diag(NullShift)*ones(3,n))';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Wideband "White" Noise %%
WhiteNoise = AddWN*(G_OutNoise*randn(n,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Gauss-Markov Correlated Noise %%
% Set up a state-space representation of noise %
a_gyro = -1/G_tau;
b_gyro = 1;
c_gyro = 1;
d_gyro = 0;

Q = 2*(G_IRb)^2/G_tau;  % Driving Noise White Power Spectral Density
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Correlate Noise %%
for k=2:n
    GyroCorrelated(k,:) = AddGM*(adg*GyroCorrelated(k-1,:)' + ug(k-1,:)')';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create Scale Factors %%
% Randomize if Monte Carlo being run, use set value for single run.
if AddSF == 0
    G_SF = [0; 0; 0];
else
    if const.MCLength > 1
        G_SF = 0.017*randn(3,1);
        G_SF
    else
        G_SF = 0.017*randn(3,1);
        G_SF
    end
end
    
%% Combine All Gyro Errors %%
Gyro_Corrupt = (diag(1+G_SF)*Gyro_Good')' + GyroCorrelated + WhiteNoise + NullShift;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end