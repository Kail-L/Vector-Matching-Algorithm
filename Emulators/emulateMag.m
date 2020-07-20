function [Mag_Corrupt, C_mag, Bias_mag] = emulateMag(t,Mag_Good,const)
%EMULATEMAG Generates a corrupted Magnetometer vector based on magnetometer
%characteristics defined in constants_struct.m.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAG_CORRUPT = emulateMag(t,Mag_Good,const) takes in a "clean" or true
% magnetometer measurement and corrupts it with potentially 4 sources of bias:
%
% 1.) Hard iron biases
% 2.) Soft iron biases
% 3.) Wide band "white" noise
% 4.) Scale factor errors
%
% The choice of biases used is dependent on a flag within the const
% structure input into the function. Note that all computation occurs in
% Gauss.
%
% SOURCES:
% D. Gebre-Egziabher, G.H. Elkaim, J.D. Powell, and B.W. Parkinson,
% "Calibration of Strapdown Magnetometers in Magnetic Field Domain,"
% ASCE Journal of Aerospace Engineering, Vol. 19, No. 2, April 2006, pp.
% 87-102. 
%
% D. Gebre-Egziabher, "Magnetometer Autocalibration Leveraging Measurement
% Locus Constraints", AIAA Journal of Aircraft, Vol. 44, No. 4, July -
% August 2007
%
% The model I am using is of the following form
%
%    M_corrupt = C*(M_Good + M_HardIron + M_WhiteNoise)
%
% where C is a 3x3 matrix of scale factor errors, missalignments, and soft
% iron biases. C can be expanded to be
%
%                       C = C_s*C_eta*C_alpha 
% where 
%
%                       C_s = [1 + sf_1 0 0
%                              0 1 + sf_2 0
%                              0 0 1 + sf_3]
%
%                       C_eta = [1 eta_3 -eta_2
%                                -eta_3 1 eta_1
%                                eta_y -eta_x 1]
%
%                       C_alpha = [1 + alpha_11 alpha_12 alpha_13
%                                  alpha_21 1 + alpha_22 alpha_23
%                                  alpha_31 alpha_32 1 + alpha_33]
% 
% C_s corresponds to a matrix of scale factor errors, C_eta corresponds to
% a matrix of small misaligments (similar to a DCM), and C_alpha
% corresponds to a matrix of induction parameters that are used to
% represent soft iron biases as the vehicle rotates in space. For the
% purposes of this code, we assume that there is no misalignement between
% the TAM and the body axes of the vehicle, and as such C_eta is identity.
%
% INPUT PARAMETERS:
% Mag_Good = nx3 array of true or "clean" mag measurements, in mGauss
% const = constants structure with relevant mag parameters
%
% OUTPUT PARAMETERS: 
% Mag_Corrupt = nx3 array of biased or "noisy" mag measurements, in mGauss
%
% VARIABLES:
% HardIron = nx3 array of hard iron constants that are added to clean mag
% WhiteNoise = nx3 array of random gaussian variables multiplied by Output
%              noise standard deviation, M_measNoise
% SoftIron = nx3 array of soft iron variables that are added to clean mag
%
% MagSF = Scale Factored mag outputs
%
% Kail Laughlin
% Updated 7/15/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define some constants %%
n = length(Mag_Good);             % Length of mag vector
dt = mean(diff(t));               % Time between mag measurements

% Pull constants from constants_struct.m %
M_OutNoise = const.M_measNoise;   % Mag sensor measurment/output noise, mGauss
M_HI = const.M_HardIron;          % Mag hard iron bias std, mGauss
M_a11 = const.M_alpha11;          % Soft iron induction parameter, 1,1 direction
M_a12 = const.M_alpha12;          % Soft iron induction parameter, 1,2 direction
M_a13 = const.M_alpha13;          % Soft iron induction parameter, 1,3 direction
M_a21 = const.M_alpha21;          % Soft iron induction parameter, 2,1 direction
M_a22 = const.M_alpha22;          % Soft iron induction parameter, 2,2 direction
M_a23 = const.M_alpha23;          % Soft iron induction parameter, 2,3 direction
M_a31 = const.M_alpha31;          % Soft iron induction parameter, 3,1 direction
M_a32 = const.M_alpha32;          % Soft iron induction parameter, 3,2 direction
M_a33 = const.M_alpha33;          % Soft iron induction parameter, 3,3 direction
M_ScaleFactor = const.M_SF;       % Mag scale factor errors, percent
AddHI = const.M_AddHardIron;      % Boolean to add hard iron biases 
AddSI = const.M_AddSoftIron;      % Boolean to add soft iron biases
AddSF = const.M_AddScaleFactor;   % Boolean to add scale factor
AddWN = const.M_AddWhiteNoise;    % Boolean to add white noise
UseSetVals = const.M_UseSetVals;  % Boolean to use set scale factors and hard iron from constants structure 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Hard Iron Bias Errors %%
if UseSetVals==1
    HardIron = AddHI*(const.M_HIVal');
    HardIron = (diag(HardIron)*ones(3,n))';
else
    HardIron = AddHI*(M_HI*randn(1,3));
    HardIron = (diag(HardIron)*ones(3,n))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Generate Wideband "White" Noise %%
WhiteNoise = AddWN*(M_OutNoise*randn(n,3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Scale Factors %%

if UseSetVals == 1
    M_ScF = const.M_SFVal';
    if AddSF ==1
        M_ScF
    end
else
    M_ScF = M_ScaleFactor*randn(3,1);
    if AddSF ==1
        M_ScF
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Add Soft Iron Biases %%
if AddSI ==1
    C_alpha = [1 + M_a11 M_a12 M_a13;
               M_a21 1 + M_a22 M_a23;
               M_a31 M_a32 1 + M_a33];
else
    C_alpha = eye(3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Combine All Gyro Errors %%

% Mag_Corrupt = (diag(1+AddSF*M_ScF)*C_alpha*Mag_Good')' + HardIron + WhiteNoise;
Mag_Corrupt = ((diag(1+AddSF*M_ScF)*C_alpha)*(Mag_Good + HardIron + WhiteNoise)')';
C_mag = (diag(1+AddSF*M_ScF)*C_alpha);
Bias_mag = HardIron;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end