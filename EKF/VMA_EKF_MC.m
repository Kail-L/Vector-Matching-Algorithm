%% Vector Matching Algorithm Extended Kalman Filter (EKF) for CubeSat %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to solve attitude estimation problem for sattelite truth data
% simulated from SatAttEst.m.
%
% Kail Laughlin
% 7/9/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set options for quadprog %
optionsQP =  optimset('Algorithm','active-set','Display','off'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Constants & SetUp Data %%

% Data from post_processing.m
%       IMU = [t, Gyro_Good, Gyro_Corrupt]
%       TrueAng = [t psi psidot theta thetadot phi phidot]
%       Mag = [B_ECI B_S]
%       Sun = [S_ECI S_S S_S_noisy]

% Constants %
constants_struct
I_b = const.I_b;                    % CubeSat Moment of Inertia
tau = const.G_tau;                  % gyro correlation time for gauss markov, s
IRb = const.G_inRunBiasSTD*d2r/3600;% In-Run bias, rad/s
ARW = sqrt(const.G_ARW^2/3600)*d2r; % rad/sqrt(sec)

% Determine which noise to populate Q_w with %
if const.G_UseARW == 0
    OutNoise = const.G_measNoise*d2r;   % rad/s
else
    OutNoise = ARW;
end

UseMeas = const.UseMeas;            % Boolean to say if we should perform measurements
Magvar = const.M_measNoise;         % Magnetometer variance, mgauss
Sunvar = const.S_measNoise*d2r;     % Sun sensor variance, rad
ASF = const.G_AddScaleFactor;       % Boolean used in F matrix to add scale factor
ESF = const.EstimateSF;             % Boolean used to estimate scale factor
SFnoise = const.G_SFNoise;          % Scale factor noise, ppm
SFLB=const.G_ScaleFactorLB;         % Lower bound scale factor matrix estimate for quadprog
SFUB=const.G_ScaleFactorUB;         % Upper bound on scale factor matrix estimate for quadprog
IFlag = const.IQ;                   % Flag to determine quaternion initialization
SR = const.MagSamp;                 % Measurement update rate
CSg = const.CSg;                    % Misalignment DCM between gyro and body 
alpha = const.alpha;                % Tuning Parameter for G matrix 
beta = const.beta;                  % Tuning Parameter for G matrix 
gamma = const.gamma;                % Tuning Parameter for G matrix

% NOTE: The tuning parameters above are used to ensure our estimate
% converges in an acceptable amount of time. However, the larger the tuning
% parameters, the larger the covariance bound will be on the estimate.
% Typically, there is a trade off between convergence speed and accuracy.

% Data %
t = IMU(:,1);                % time vector, s
drl=length(t);               % Looping variable for integration
dt=mean(diff(t));            % Time step of gyro timing data

uECI=Mag(:,1:3)*T2mg;        % Clean Mag field vector in inertial frame, mGauss
Mag_Clean_S=Mag(:,4:6)*T2mg; % Clean Mag field vector in body frame, mGauss        

sECI=Sun(:,1:3);             % Unit sun vector measurement in inertial frame
Sun_Clean_S=Sun(:,4:6);      % Unit sun vector measurement in body frame
sS=Sun(:,7:9);               % Unit sun vector measurment in body frame

Om_SE_good=IMU(:,2:4);       % Truth angular rate values 
psi_S_ECI=TrueAng(:,2)*r2d;  % Truth value of psi, deg
the_S_ECI=TrueAng(:,4)*r2d;  % Truth value of theta, deg
phi_S_ECI=TrueAng(:,6)*r2d;  % Truth value of phi, deg
q_S_ECI=q_S_ECI_true;        % True quaternion (note, this is unecessary but placed here for clarity)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Monte Carlo Storage %
qhat0_S_ECI_MC = zeros(drl,const.MCLength);
qhat1_S_ECI_MC = zeros(drl,const.MCLength);
qhat2_S_ECI_MC = zeros(drl,const.MCLength);
qhat3_S_ECI_MC = zeros(drl,const.MCLength);
q_SS_hat1_MC = zeros(drl,const.MCLength);
q_SS_hat2_MC = zeros(drl,const.MCLength);
q_SS_hat3_MC = zeros(drl,const.MCLength);
std_qSS1_MC = zeros(drl,const.MCLength);
std_qSS2_MC = zeros(drl,const.MCLength);
std_qSS3_MC = zeros(drl,const.MCLength);

b1_MC = zeros(drl,const.MCLength);
b2_MC = zeros(drl,const.MCLength);
b3_MC = zeros(drl,const.MCLength);
del_b1_MC = zeros(drl,const.MCLength);
del_b2_MC = zeros(drl,const.MCLength);
del_b3_MC = zeros(drl,const.MCLength);
std_db1_MC = zeros(drl,const.MCLength);
std_db2_MC = zeros(drl,const.MCLength);
std_db3_MC = zeros(drl,const.MCLength);

sf1_MC = zeros(drl,const.MCLength);
sf2_MC = zeros(drl,const.MCLength);
sf3_MC = zeros(drl,const.MCLength);
del_sf1_MC = zeros(drl,const.MCLength);
del_sf2_MC = zeros(drl,const.MCLength);
del_sf3_MC = zeros(drl,const.MCLength);
std_dsf1_MC = zeros(drl,const.MCLength);
std_dsf2_MC = zeros(drl,const.MCLength);
std_dsf3_MC = zeros(drl,const.MCLength);

phihat_S_ECI_MC = zeros(drl,const.MCLength);
thehat_S_ECI_MC = zeros(drl,const.MCLength);
psihat_S_ECI_MC = zeros(drl,const.MCLength);
psi_SS_hat_MC = zeros(drl,const.MCLength);
the_SS_hat_MC = zeros(drl,const.MCLength);
phi_SS_hat_MC = zeros(drl,const.MCLength);
std_psiSS_MC = zeros(drl,const.MCLength);
std_theSS_MC = zeros(drl,const.MCLength);
std_phiSS_MC = zeros(drl,const.MCLength);

Epsilon_MC = zeros(drl,const.MCLength);
std_epsilon_MC = zeros(drl,const.MCLength);

Obs_G_rank_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv_MC = zeros(length(t_mag),const.MCLength);
Obs_G_cond_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv1_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv2_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv3_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv4_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv5_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv6_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv7_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv8_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv9_MC = zeros(length(t_mag),const.MCLength);
L_sv_MC = zeros(length(t_mag),const.MCLength);
Q_sv_MC = zeros(length(t_mag),const.MCLength);
Obs_G_rank_true_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv_true_MC = zeros(length(t_mag),const.MCLength);
Obs_G_cond_true_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv1_true_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv2_true_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv3_true_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv4_true_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv5_true_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv6_true_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv7_true_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv8_true_MC = zeros(length(t_mag),const.MCLength);
Obs_G_sv9_true_MC = zeros(length(t_mag),const.MCLength);
L_sv_true_MC = zeros(length(t_mag),const.MCLength);
Q_sv_true_MC = zeros(length(t_mag),const.MCLength);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Begin Monte Carlo %%

tic
for MC=1:const.MCLength
%% Set up Data and initialize matrices for  s p e e d %%
% Initialize Matrices %
qhat_S_ECI = zeros(drl,4);         % Estimate quaternion
q_SS_hat = zeros(drl,4);           % Quaternion errors from body frame to estimate body frame
psihat_S_ECI=zeros(drl,1);         % Estimate value of psi, deg
thehat_S_ECI=zeros(drl,1);         % Estimate value of theta, deg
phihat_S_ECI=zeros(drl,1);         % Estimate value of phi, deg
psi_SS_hat = zeros(drl,1);         % Error in psi, deg
the_SS_hat = zeros(drl,1);         % Error in theta, deg
phi_SS_hat = zeros(drl,1);         % Error in phi, deg

del_q_SS = zeros(drl,3);           % Quaternion vector estimation errors
del_b = zeros(drl,3);              % Bias estimation errors
del_sf = zeros(drl,3);             % Scale factor estimation errrors

b = zeros(drl,3);                  % Storage vector of biases
b_out = [0 0 0]'*d2r;              % Initial bias estimate on gyros
b(1,:) = b_out';                   % Initial bias estimate
sf = zeros(drl,3);                % Storage vector of scale factors
sf_out=[0 0 0]';                  % Initial inverse SF estimate on gyros
sf(1,:) = sf_out';                  % Initial inverse scale factor estimate

Epsilon = zeros(drl,1);            % Angular error between estimate and true body 3 axes
uS_rec = zeros(length(t_mag),3);   % Magnetic field in body frame measured                              
uS_comp = zeros(length(t_mag),3);  % Magnetic field in body frame computed with true attitude
sS_rec = zeros(length(t_sun),3);   % Sun vector in body frame measured
sS_comp = zeros(length(t_sun),3);  % Sun vector in body frame computed with true attitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% State Covariance and Storage %
if ESF ==1
    P_stor = zeros(drl,9);         % Covariance Storage vector
else
    P_stor = zeros(drl,6);
end
P_eul = zeros(4,4);                % Covariance matrix used to convert quat covariance to eul covariance

std_psiSS = zeros(drl,1);
std_theSS = zeros(drl,1);
std_phiSS = zeros(drl,1);

std_qSS1 = zeros(drl,1);
std_qSS2 = zeros(drl,1);
std_qSS3 = zeros(drl,1);

std_db1 = zeros(drl,1);
std_db2 = zeros(drl,1);
std_db3 = zeros(drl,1);

std_dsf1 = zeros(drl,1);
std_dsf2 = zeros(drl,1);
std_dsf3 = zeros(drl,1);

std_epsilon = zeros(drl,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Observability Values for VMA_SO_Test.m %
if ESF ==1
    Fd_TM = eye(9);                 % Initial State transition matrix of F
else
    Fd_TM = eye(6);
end

Obs = [];                           % Initial observability matrix
Obs_rank = zeros(length(t_mag),1);  % Storage vector for rank of Observability matrix
Obs_G_rank = zeros(length(t_mag),1);% Observability grammian rank check

H_Test = cell(length(t_mag),1);       % Storage cell of measurement matrices
H_Test_true = cell(length(t_mag),1);
H_row = 0;                          % Initial value of number of rows of H
H_col = 0;                          % Initial value of number of columns of H
PHI_Test = cell(drl,1);               % Storage cell for current state transition matrix
PHI_Test_true = cell(drl,1);               % Storage cell for current state transition matrix
PHI_Test{1} = Fd_TM;                  % Initial state transistion matrix
PHI_Test_true{1} = Fd_TM;                  % Initial state transistion matrix
PhiSize = size(Fd_TM);              % Size of state transition matrix
PHI_row = PhiSize(1);               % Rows of state transition matrix
PHI_col = PhiSize(2);               % Columns of state transition matrix
Rz_Test = cell(drl,1);                % Storage cell for measurement covariance matrix
Rz_Test_true = cell(drl,1);                % Storage cell for measurement covariance matrix
Qd_Test = cell(drl,1);                % Storage cell for process noise covariance matrix
Qd_Test_true = cell(drl,1);                % Storage cell for process noise covariance matrix
Qd_row = 0;                         % Initial rows of process noise covariance matrix
Qd_col = 0;                         % Initial columns of process noise covariance matrix
TU_num = drl;                       % Number of time updates
TU_time = t;                        % Time for time updates
MU_num = length(t_mag);             % Number of measurement updates
MU_time = t_mag;                    % Time for measurement updates
tu_per_mu = const.GyroSamp/const.MagSamp;% Ratio of time updates to measurement updates
state_cov_his = cell(drl,1);        % Storage cell for state covariance matrix P

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% EKF Initialization %%

disp(['Monte Carlo Run: ',num2str(MC),'/',num2str(const.MCLength)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corrupt Gyro Parameters %
[Om_SE_noisy,G_SF,B_true] = emulateGyro(t,Om_SE_good,const);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corrupt Mag Parameters %
if const.UseSGPB == 1
    uS = B_S*T2mg;      % Convert to mGauss
else
    uS = emulateMag(t_mag,Mag_Clean_S,const);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial State Covariance %
P0=const.P0;
if ESF ==1 % If we are estimating scale factor
    P_stor(1,:) = [P0(1,1) P0(2,2) P0(3,3) P0(4,4) P0(5,5) P0(6,6) P0(7,7) P0(8,8) P0(9,9)];
    std_qSS1(1) = sqrt(P0(1,1));
    std_qSS2(1) = sqrt(P0(2,2));
    std_qSS3(1) = sqrt(P0(3,3));
    std_db1(1) = sqrt(P0(4,4));
    std_db2(1) = sqrt(P0(5,5));
    std_db3(1) = sqrt(P0(6,6));
    std_dsf1 = sqrt(P0(7,7));
    std_dsf2 = sqrt(P0(8,8));
    std_dsf3 = sqrt(P0(9,9));
else % If we are NOT estimating scale factor
    P_stor(1,:) = [P0(1,1) P0(2,2) P0(3,3) P0(4,4) P0(5,5) P0(6,6)];        
    std_qSS1(1) = sqrt(P0(1,1));
    std_qSS2(1) = sqrt(P0(2,2));
    std_qSS3(1) = sqrt(P0(3,3));
    std_db1(1) = sqrt(P0(4,4));
    std_db2(1) = sqrt(P0(5,5));
    std_db3(1) = sqrt(P0(6,6));
end

state_cov_his{1}=P0;             % Store for SO test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensor Covariances %
Rmag = ((Magvar)^2)*eye(3);    % Magnetometer covariance matrix, mGauss
Rsun = ((Sunvar)^2)*eye(3);    % Sun Sensor covariance matrix
if const.UseSun==1
    Rz(1:3,1:3) = Rmag;        % State Measurement covariance matrix
    Rz(4:6,4:6) = Rsun;
else
    Rz = Rmag;                 % State Measurement covariance matrix
end

Rz_Test{1} = Rz;       % Store for SO test, convert to mGauss
Rz_Test_true{1} = Rz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process Noise Covariance %
if ESF ==1
    Qw = [(OutNoise)^2*eye(3) zeros(3,3) zeros(3,3);   % Gyro Output noise
          zeros(3,3) 2*(IRb)^2/tau*eye(3) zeros(3,3);  % Gyro in-run bias repeat.
          zeros(3,3) zeros(3,3) (SFnoise/1e6)^2*eye(3)];    % Scale factor noise (Based on gyro params sent by Vibhor)
else
    Qw = [(OutNoise)^2*eye(3) zeros(3,3);   % Gyro Output noise
          zeros(3,3) 2*(IRb)^2/tau*eye(3)]; % Gyro in-run bias repeat.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Quaternion %%
q_S_ECI_0 = q_S_ECI(1,:);        % Initial truth quaternion

if IFlag == 0                 % IFlag == 0: Use Truth quat as initial estimate
    qhat_S_ECI_0 = q_S_ECI_0';  
elseif IFlag == 1             % IFlag == 1: Use Max angular error allowed, deg
    MaxAngErr=const.ME;       % Maximum angular error allowed, deg
    ErrEul = [MaxAngErr*(2*rand-1);
              MaxAngErr*(2*rand-1);
              MaxAngErr*(2*rand-1)];
    ErrEul = ErrEul*d2r;
    C_SECIErr=eul2dcm(ErrEul);
    Chat_S_ECI_0 = C_SECIErr*Cba_0;
    qhat_S_ECI_0 = dcm2quat(Chat_S_ECI_0)';
elseif IFlag == 2             % IFlag == 2: Use completely random quat
    q0=rand(1,1);             % Random quaternion
    q1=2*rand(1,1) - 1;
    q2=2*rand(1,1) - 1;
    q3=2*rand(1,1) - 1;
    qhat_S_ECI_0 = [q0 q1 q2 q3]';
    qhat_S_ECI_0 = qhat_S_ECI_0/quatnorm(qhat_S_ECI_0);
elseif IFlag == 3             % IFlag == 3: Use set of error Euler angles
    ErrEul=const.ErrEul;      % Euler angle error
    C_SECIErr=eul2dcm(ErrEul);
    Chat_S_ECI_0 = C_SECIErr*Cba_0;
    qhat_S_ECI_0 = dcm2quat(Chat_S_ECI_0)';
end

qhat_S_ECI(1,:) = qhat_S_ECI_0';    % Initial quaternion estimate based on above
if qhat_S_ECI(1,1)<0                % Enforce positive scalar
    qhat_S_ECI(1,:)=-qhat_S_ECI(1,:);
end

% Initial Euler angle estimate based on above quaternion %
EE_hat_0 = dcm2eul(quat2dcm(qhat_S_ECI_0'))*180/pi;
psihat_S_ECI_0 = EE_hat_0(1); % yaw
thehat_S_ECI_0 = EE_hat_0(2); % pitch
phihat_S_ECI_0 = EE_hat_0(3); % roll

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Initial Angular Error, Epsilon %%
Chat_S_ECI_0 = quat2dcm(qhat_S_ECI_0');     
C_SS_hat_0 = quat2dcm(q_S_ECI_0)*inv(Chat_S_ECI_0);
vS = [0 0 1]';          % 3-axis of S frame
vS_hat = C_SS_hat_0*vS; % 3-axis of S_hat frame
epsilon = acos(C_SS_hat_0(3,3))*180/pi;

disp(['Initial Epsilon (Angular Error Between S3 and S3_hat): ',num2str(epsilon),' degrees'])
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize quaternion error, body frame %%
q_SS_hat_0 = quatmult(q_S_ECI_0',quatInverse(qhat_S_ECI_0'));
q_SS_hat_0 = q_SS_hat_0/norm(q_SS_hat_0);
q_SS_hat(1,:) = q_SS_hat_0';
del_q_SS(1,:) = q_SS_hat(1,2:4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Yaw, Pitch, Roll, degrees %
eulhat_S_ECI_0 = dcm2eul(Chat_S_ECI_0);
psihat_S_ECI(1) = eulhat_S_ECI_0(1)*r2d;
thehat_S_ECI(1) = eulhat_S_ECI_0(2)*r2d;
phihat_S_ECI(1) = eulhat_S_ECI_0(3)*r2d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run EKF %%
mu_counter=2;           % Measurement update counter 
disp('Running EKF...')
fprintf('\n')
for k = 2:drl
%%%%%%%% Time Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ESF ==1
        %Estimated Trajectory%
        iSFmat = inv(diag([1+sf_out(1) 1+sf_out(2) 1+sf_out(3)]));
        Om_hat_SE = iSFmat*CSg*(Om_SE_noisy(k,:)'-b_out);% Estimate angular velocity, body frame

        F = [-sk(Om_hat_SE) -1/2*iSFmat*CSg -1/2*diag(CSg*(Om_SE_noisy(k,:)'-b_out)); 
            zeros(3) (-1/tau)*eye(3) zeros(3);
            zeros(3) zeros(3) zeros(3)];

        G = [-1/2*alpha*iSFmat*CSg zeros(3) zeros(3);         
            zeros(3) beta*CSg zeros(3);
            zeros(3) zeros(3) gamma*eye(3)];
     
        %Truth Trajectory%
        iSFmat_true = inv(diag([1+G_SF(1) 1+G_SF(2) 1+G_SF(3)]));
        Om_SE = Om_SE_good(k,:)';
        Ftrue = [-sk(Om_SE) -1/2*(iSFmat_true)*CSg -1/2*diag(CSg*(Om_SE_noisy(k,:)'-B_true(k,:)')); 
                 zeros(3) (-1/tau)*eye(3) zeros(3);
                 zeros(3) zeros(3) zeros(3)];   
        Gtrue = [-1/2*iSFmat_true*alpha*CSg zeros(3) zeros(3);         
                 zeros(3) beta*CSg zeros(3);
                 zeros(3) zeros(3) gamma*eye(3)];
    else
        %Estimated Trajectory%
        Om_hat_SE = Om_SE_noisy(k,:)'-b_out;             % Estimate angular velocity, body frame
        F = [-sk(Om_hat_SE) -1/2*eye(3); 
            zeros(3) (-1/tau)*eye(3)];

        G = [-1/2*alpha*CSg zeros(3);         
            zeros(3) beta*eye(3)];
       
        %Truth Trajectory%
        Om_SE = Om_SE_good(k,:)';
        Ftrue = [-sk(Om_SE) -1/2*eye(3);
                  zeros(3) (-1/tau)*eye(3)];
        Gtrue = [-1/2*alpha*CSg zeros(3);
                  zeros(3) beta*eye(3)];
    end
    
    FullRank = max(size(F));

    % Propagate deltas during no measurement
    del_b(k,:)=del_b(k-1,:);
    del_sf(k,:)=del_sf(k-1,:);
    del_q_SS(k,:)=del_q_SS(k-1,:);

    % Covariance Update 
    Phi = expm(F*dt);           % Discreet time F matrix (Phi matrix)
    PhiTrue = expm(Ftrue*dt);   % Discreet time Ftrue matrix
    Cd = disrw(F,G,dt,Qw);      % Discreet time noise mapping matrix (Gamma matrix)
    CdTrue = disrw(Ftrue,Gtrue,dt,Qw);% Discreet time noise mapping matrix
    P = Phi*P0*Phi' + Cd;       % Updated covariance
    PhiSize = size(Phi);
    CdSize = size(Cd);        

    % Estimate Quaternion Update %
    qhat_S_ECI(k,:) = quatIntegrate((Om_hat_SE)',qhat_S_ECI(k-1,:)',dt);
    qhat_S_ECI(k,:) = qhat_S_ECI(k,:)/quatnorm(qhat_S_ECI(k,:));
    if (qhat_S_ECI(k,1) < 0)         % Enforce Positive Scalar 
        qhat_S_ECI(k,:) = -qhat_S_ECI(k,:);
    end

    % Enforce Positive Scalar on True Quaternion %
    if q_S_ECI(k,1) < 0              
        q_S_ECI(k,:) = -q_S_ECI(k,:);
    end

    % Propagate bias, sf %
    b(k,:)=b(k-1,:);
    sf(k,:)=sf(k-1,:); 

    % Store Variables for SO test %
    PHI_Test{k} = Phi;
    PHI_Test_true{k} = PhiTrue;
    PHI_row = PhiSize(1);
    PHI_col = PhiSize(2);
    Qd_Test{k} = Cd;
    Qd_Test_true{k} = CdTrue;
    Qd_row = CdSize(1);
    Qd_col = CdSize(2);

%%%%%%%% Measurement Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if UseMeas && mu_counter<=length(t_mag)&&t_mag(mu_counter)==t(k)

        % Current Estimate DCM %
        Chat_S_ECI = quat2dcm(qhat_S_ECI(k,:));
        C_S_ECI = quat2dcm(q_S_ECI(k,:));
        
        % Non-Normalized Mag Field model & measurement 
        US = uS(mu_counter,:)';     uS_rec(mu_counter,:) = US';
        UECI = uECI(mu_counter,:)';

        % Normalized Sun Sensor model & measurement (if used) %

        SS = sS(mu_counter,:)';     sS_rec(mu_counter,:) = SS';
        SECI = sECI(mu_counter,:)';  

        % Measurement Matrix %
        if const.UseSun==1
            H = [-2*sk(Chat_S_ECI*UECI) zeros(3,FullRank-3);
                 -2*sk(Chat_S_ECI*SECI) zeros(3,FullRank-3)]; 
            Htrue = [-2*sk(C_S_ECI*UECI) zeros(3,FullRank-3);
                     -2*sk(C_S_ECI*SECI) zeros(3,FullRank-3)];
        else
            H = [-2*sk(Chat_S_ECI*UECI) zeros(3,FullRank-3)];
            Htrue = [-2*sk(C_S_ECI*UECI) zeros(3,FullRank-3)];
        end

        % Store values for stochastic observability test %
        H_Test{mu_counter} = H;   
        H_Test_true{mu_counter} = Htrue;
        HSize = size(H);
        H_row = HSize(1);
        H_col = HSize(2);

        mag_S = Chat_S_ECI*UECI;
        sun_S = Chat_S_ECI*SECI;
        dmag = mag_S-US;
        dsun = sun_S-SS;
        if const.UseSun==1
            dz = [dmag;
                  dsun];
        else 
            dz = dmag;
        end

%%%%%%%%%%% Run EKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        Rz_Test{mu_counter} = Rz; 
        Rz_Test_true{mu_counter} = Rz;
        SizeP = size(P);
        K = (P*H')/(H*P*H' + Rz);       % Kalman Filter Gain 
        P = (eye(SizeP(1)) - K*H)*P;    % Covariance update 
        del_x = K*dz;                   % State Error Vector 

%%%%%%%%%%% Constrain Estimates if scale factor used %%%%%%%%%%%%%%%%%%%%%%            
        if ESF ==1
            if const.constrainEst == 1
            % Constrain using QuadProg
                D=[];
                d=[];
                Aeq  = [0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0;
                        0 0 0 0 0 0 0 0 0];
                beq = [0 0 0 0 0 0 0 0 0];                   
                lb_con=[-inf -inf -inf -inf -inf -inf 1+SFLB 1+SFLB 1+SFLB]'-[0 0 0 0 0 0 (1+sf_out')]';
                ub_con=[inf inf inf inf inf inf 1+SFUB 1+SFUB 1+SFUB]'-[0 0 0 0 0 0 (1+sf_out')]';
                x0 = del_x;
                Hqp = inv(P);
                Hqp = (Hqp+Hqp')./2;
                xdbar=del_x';
                ft=-xdbar*Hqp;
                f = ft';
                xbar=quadprog(Hqp,f,D,d,Aeq,beq,lb_con,ub_con,x0,optionsQP);
            else
                xbar = del_x;
            end
        end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Store Error State %
        if ESF == 1
            del_q_SS(k,:) = xbar(1:3);
            del_b(k,:) = xbar(4:6);
            del_sf(k,:) = xbar(7:9);
        else
            del_q_SS(k,:)= del_x(1:3)';
            del_b(k,:)= del_x(4:6)';   
        end

        % Update States %
        q_SS_hat(k,:)=[1 del_q_SS(k,:)];
        b_out = b_out + del_b(k,:)';
        b(k,:) = b_out;

        sf_out = sf_out + del_sf(k,:)';
        sf(k,:) = sf_out';
        qhat_S_ECI(k,:) = quatmult(q_SS_hat(k,:)',qhat_S_ECI(k,:)');
        qhat_S_ECI(k,:) = qhat_S_ECI(k,:)'/quatnorm(qhat_S_ECI(k,:));

        % Update Measurement Counter %
        mu_counter = mu_counter+1;
    end

%%%%%%%% Package Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    % Enforce Positive Scalar %
    if (qhat_S_ECI(k,1) < 0)
        qhat_S_ECI(k,:) = -qhat_S_ECI(k,:);
    end

    % Error Quaternion Update%
    if const.UseSGPB == 1
        q_SS_hat(k,:) = quatmult(q_S_ECI(k,:)',quatInverse(qhat_S_ECI(k,:)));
    else
        q_SS_hat(k,:) = quatmult(q_S_ECI(k,:)',quatInverse(qhat_S_ECI(k,:)));
    end
    q_SS_hat(k,:) = q_SS_hat(k,:)/quatnorm(q_SS_hat(k,:));

    if (q_SS_hat(k,1) < 0)         % Enforce Positive Scalar 
        q_SS_hat(k,:) = -q_SS_hat(k,:);
    end


    % Store state covariance history
    state_cov_his{k}=P;

    % Store Covariance Estimate %
    if ESF ==1
        P_stor(k,:) = [P(1,1) P(2,2) P(3,3) P(4,4) P(5,5) P(6,6) P(7,7) P(8,8) P(9,9)];
        P0 = P;
    else
        P_stor(k,:) = [P(1,1) P(2,2) P(3,3) P(4,4) P(5,5) P(6,6)];
        P0 = P;           
    end

    % Updated DCM from measurement quaternion
    Chat_S_ECI = quat2dcm(qhat_S_ECI(k,:));  

    % Error DCM %
    C_SS_hat = quat2dcm(q_SS_hat(k,:));

    % Angular Error %
    vS_hat = C_SS_hat*vS;
    Epsilon(k) = acos(C_SS_hat(3,3))*180/pi;

    % Standard deviation of quaternions and biases %
    std_qSS1(k) = sqrt(P(1,1));
    std_qSS2(k) = sqrt(P(2,2));
    std_qSS3(k) = sqrt(P(3,3));
    std_db1(k) = sqrt(P(4,4));
    std_db2(k) = sqrt(P(5,5));
    std_db3(k) = sqrt(P(6,6));
    if ESF ==1
        std_dsf1(k) = sqrt(P(7,7));
        std_dsf2(k) = sqrt(P(8,8));
        std_dsf3(k) = sqrt(P(9,9));
    end

    % Estimated Euler Angles, deg %
    eul_hat_ba = dcm2eul(Chat_S_ECI);
    psihat_S_ECI(k,1) = eul_hat_ba(1)*r2d;
    thehat_S_ECI(k,1) = eul_hat_ba(2)*r2d;
    phihat_S_ECI(k,1) = eul_hat_ba(3)*r2d;

    % Estimated Euler Angle Errors %
    err_eul_hat = dcm2eul(C_SS_hat);
    psi_SS_hat(k) = err_eul_hat(1)*r2d;
    the_SS_hat(k) = err_eul_hat(2)*r2d;
    phi_SS_hat(k) = err_eul_hat(3)*r2d;

    % Standard deviation of euler angles %
    Gq=QuatJoc(q_SS_hat(k,:));
    P_eul=P(1:3,1:3); 
    Peul=Gq(:,2:4)*P_eul*Gq(:,2:4)';
    std_psiSS(k) = sqrt(Peul(3,3))*r2d;
    std_theSS(k) = sqrt(Peul(2,2))*r2d;
    std_phiSS(k) = sqrt(Peul(1,1))*r2d;

    % Standard deviation of epsilon %
    Eq = EpsJoc(q_SS_hat(k,:));
    P_eps=P(1:3,1:3);
    Peps=Eq(2:4)*P_eps*Eq(2:4)';
    std_epsilon(k) = sqrt(Peps)*r2d;
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Stochastic Observability Test %%
if const.RunSOTest ==1
    % Estimate Trajectory %
    PHI_SO = PHI_Test;
    Qd_SO = Qd_Test;
    H_SO = H_Test;
    Rz_SO = Rz_Test;
    VMA_SO_Test                         % Run SO test
    Obs_G_rank_MC(:,MC) = r_OBS_gram';           % Rank of observability gramian
    Obs_G_cond_MC(:,MC) = c_OBS_gram';           % Condition # of observability gramian
    Obs_G_sv1_MC(:,MC) = OBS_gram_sv(1,:)';      % 1st singular value of obsv. gram.
    Obs_G_sv2_MC(:,MC) = OBS_gram_sv(2,:)';      % 2nd singular value of obsv. gram.
    Obs_G_sv3_MC(:,MC) = OBS_gram_sv(3,:)';      % 3rd singular value of obsv. gram.
    Obs_G_sv4_MC(:,MC) = OBS_gram_sv(4,:)';      % 4th singular value of obsv. gram.
    Obs_G_sv5_MC(:,MC) = OBS_gram_sv(4,:)';      % 5th singular value of obsv. gram.
    Obs_G_sv6_MC(:,MC) = OBS_gram_sv(6,:)';      % 6th singular value of obsv. gram.
    if ESF==1
        Obs_G_sv7_MC(:,MC) = OBS_gram_sv(7,:)';  % 7th singular value of obsv. gram.
        Obs_G_sv8_MC(:,MC) = OBS_gram_sv(7,:)';  % 8th singular value of obsv. gram.
        Obs_G_sv9_MC(:,MC) = OBS_gram_sv(7,:)';  % 9th singular value of obsv. gram.
    end
    L_sv_MC(:,MC) = (Lambda_sv(1,:))';      % Max singular value of Lambda metric (see Chapter 4. of Kail's Thesis)
    Q_sv_MC(:,MC) = (Q_bar_sv(1,:))';        % Max singular value of Q bar metric (see Chapter 4. of Kail's Thesis)

    % Truth Trajectory %
    PHI_SO = PHI_Test_true;
    Qd_SO = Qd_Test_true;
    H_SO = H_Test_true;
    Rz_SO = Rz_Test_true;
    VMA_SO_Test                              % Run SO test
    Obs_G_rank_true_MC(:,MC) = r_OBS_gram';           % Rank of observability gramian
    Obs_G_cond_true_MC(:,MC) = c_OBS_gram';           % Condition # of observability gramian
    Obs_G_sv1_true_MC(:,MC) = OBS_gram_sv(1,:)';      % 1st singular value of obsv. gram.
    Obs_G_sv2_true_MC(:,MC) = OBS_gram_sv(2,:)';      % 2nd singular value of obsv. gram.
    Obs_G_sv3_true_MC(:,MC) = OBS_gram_sv(3,:)';      % 3rd singular value of obsv. gram.
    Obs_G_sv4_true_MC(:,MC) = OBS_gram_sv(4,:)';      % 4th singular value of obsv. gram.
    Obs_G_sv5_true_MC(:,MC) = OBS_gram_sv(4,:)';      % 5th singular value of obsv. gram.
    Obs_G_sv6_true_MC(:,MC) = OBS_gram_sv(6,:)';      % 6th singular value of obsv. gram.
    if ESF==1
        Obs_G_sv7_true_MC(:,MC) = OBS_gram_sv(7,:)';  % 7th singular value of obsv. gram.
        Obs_G_sv8_true_MC(:,MC) = OBS_gram_sv(7,:)';  % 8th singular value of obsv. gram.
        Obs_G_sv9_true_MC(:,MC) = OBS_gram_sv(7,:)';  % 9th singular value of obsv. gram.
    end
    L_sv_true_MC(:,MC) = (Lambda_sv(1,:))';       % Max singular value of Lambda metric (see Chapter 4. of Kail's Thesis)
    Q_sv_true_MC(:,MC) = (Q_bar_sv(1,:))';        % Max singular value of Q bar metric (see Chapter 4. of Kail's Thesis)
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Monte Carlo Storage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

q_SS_hat1_MC(:,MC) =  q_SS_hat(:,2);
q_SS_hat2_MC(:,MC) =  q_SS_hat(:,3);
q_SS_hat3_MC(:,MC) =  q_SS_hat(:,4);

std_qSS1_MC(:,MC) = sqrt(P_stor(:,1));
std_qSS2_MC(:,MC) = sqrt(P_stor(:,2));
std_qSS3_MC(:,MC) = sqrt(P_stor(:,3));

del_b1_MC(:,MC) = del_b(:,1);
del_b2_MC(:,MC) = del_b(:,2);
del_b3_MC(:,MC) = del_b(:,3);
std_db1_MC(:,MC) = sqrt(P_stor(:,4));
std_db2_MC(:,MC) = sqrt(P_stor(:,5));
std_db3_MC(:,MC) = sqrt(P_stor(:,6));

if ESF==1
    del_sf1_MC(:,MC) = del_sf(:,1);
    del_sf2_MC(:,MC) = del_sf(:,2);
    del_sf3_MC(:,MC) = del_sf(:,3);
    std_dsf1_MC(:,MC) = sqrt(P_stor(:,7));
    std_dsf2_MC(:,MC) = sqrt(P_stor(:,8));
    std_dsf3_MC(:,MC) = sqrt(P_stor(:,9));
    sf1_MC(:,MC) = sf(:,1);
    sf2_MC(:,MC) = sf(:,2);
    sf3_MC(:,MC) = sf(:,3);
end

qhat0_S_ECI_MC(:,MC) = qhat_S_ECI(:,1);
qhat1_S_ECI_MC(:,MC) = qhat_S_ECI(:,2);
qhat2_S_ECI_MC(:,MC) = qhat_S_ECI(:,3);
qhat3_S_ECI_MC(:,MC) = qhat_S_ECI(:,4);

b1_MC(:,MC) = b(:,1);
b2_MC(:,MC) = b(:,2);
b3_MC(:,MC) = b(:,3);
    
phihat_S_ECI_MC(:,MC) = phihat_S_ECI;
thehat_S_ECI_MC(:,MC) = thehat_S_ECI;
psihat_S_ECI_MC(:,MC) = psihat_S_ECI;

psi_SS_hat_MC(:,MC) = psi_SS_hat;  
std_psiSS_MC(:,MC) = std_psiSS;
the_SS_hat_MC(:,MC) = the_SS_hat;
std_theSS_MC(:,MC) = std_theSS;
phi_SS_hat_MC(:,MC) = phi_SS_hat;
std_phiSS_MC(:,MC) = std_phiSS;

Epsilon_MC(:,MC) = Epsilon;
std_epsilon_MC(:,MC) = std_epsilon;

% Approximate MC Runtime %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TTT=toc;
if const.MCLength>1&&MC==1
    c=clock;
    t_sim_left = datevec(datetime(c)+seconds(TTT*const.MCLength));
    d=datetime(t_sim_left);
    disp(['Approximate Monte Carlo Completion: ',datestr(d)])
    fprintf('\n')
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compute Mean of Monte Carlo Values %%
if const.plotMean ==1
    % Estimate Trajectory %
    Obs_G_rank_Mean = sum(Obs_G_rank_MC,2)./MC;      % Rank of observability gramian
    Obs_G_cond_Mean = sum(Obs_G_cond_MC,2)./MC;      % Condition # of observability gramian
    Obs_G_sv1_Mean = sum(Obs_G_sv1_MC,2)./MC;      % 1st singular value of obsv. gram.
    Obs_G_sv2_Mean = sum(Obs_G_sv2_MC,2)./MC;      % 2nd singular value of obsv. gram.
    Obs_G_sv3_Mean = sum(Obs_G_sv3_MC,2)./MC;      % 3rd singular value of obsv. gram.
    Obs_G_sv4_Mean = sum(Obs_G_sv4_MC,2)./MC;      % 4th singular value of obsv. gram.
    Obs_G_sv5_Mean = sum(Obs_G_sv5_MC,2)./MC;      % 5th singular value of obsv. gram.
    Obs_G_sv6_Mean = sum(Obs_G_sv6_MC,2)./MC;      % 6th singular value of obsv. gram.
    if ESF==1
        Obs_G_sv7_Mean = sum(Obs_G_sv7_MC,2)./MC;  % 7th singular value of obsv. gram.
        Obs_G_sv8_Mean = sum(Obs_G_sv8_MC,2)./MC;  % 8th singular value of obsv. gram.
        Obs_G_sv9_Mean = sum(Obs_G_sv9_MC,2)./MC;  % 9th singular value of obsv. gram.
    end
    L_sv_Mean = sum(L_sv_MC,2)./MC;
    Q_sv_Mean = sum(Q_sv_MC,2)./MC;

    % Truth Trajectory %
    Obs_G_rank_true_Mean = sum(Obs_G_rank_true_MC,2)./MC;      % Rank of observability gramian
    Obs_G_cond_true_Mean = sum(Obs_G_cond_true_MC,2)./MC;      % Condition # of observability gramian
    Obs_G_sv1_true_Mean = sum(Obs_G_sv1_true_MC,2)./MC;      % 1st singular value of obsv. gram.
    Obs_G_sv2_true_Mean = sum(Obs_G_sv2_true_MC,2)./MC;      % 2nd singular value of obsv. gram.
    Obs_G_sv3_true_Mean = sum(Obs_G_sv3_true_MC,2)./MC;      % 3rd singular value of obsv. gram.
    Obs_G_sv4_true_Mean = sum(Obs_G_sv4_true_MC,2)./MC;      % 4th singular value of obsv. gram.
    Obs_G_sv5_true_Mean = sum(Obs_G_sv5_true_MC,2)./MC;      % 5th singular value of obsv. gram.
    Obs_G_sv6_true_Mean = sum(Obs_G_sv6_true_MC,2)./MC;      % 6th singular value of obsv. gram.
    if ESF==1
        Obs_G_sv7_true_Mean = sum(Obs_G_sv7_true_MC,2)./MC;  % 7th singular value of obsv. gram.
        Obs_G_sv8_true_Mean = sum(Obs_G_sv8_true_MC,2)./MC;  % 8th singular value of obsv. gram.
        Obs_G_sv9_true_Mean = sum(Obs_G_sv9_true_MC,2)./MC;  % 9th singular value of obsv. gram.
    end
    L_sv_true_Mean = sum(L_sv_true_MC,2)./MC;
    Q_sv_true_Mean = sum(Q_sv_true_MC,2)./MC;

    q_SS_hat1_Mean =  sum(q_SS_hat1_MC,2)./MC;
    q_SS_hat2_Mean =  sum(q_SS_hat2_MC,2)./MC;
    q_SS_hat3_Mean =  sum(q_SS_hat3_MC,2)./MC;

    std_qSS1_Mean = sum(std_qSS1_MC,2)./MC;
    std_qSS2_Mean = sum(std_qSS2_MC,2)./MC;
    std_qSS3_Mean = sum(std_qSS3_MC,2)./MC;

    del_b1_Mean = sum(del_b1_MC,2)./MC;
    del_b2_Mean = sum(del_b2_MC,2)./MC;
    del_b3_Mean = sum(del_b3_MC,2)./MC;
    std_db1_Mean = sum(std_db1_MC,2)./MC;
    std_db2_Mean = sum(std_db2_MC,2)./MC;
    std_db3_Mean = sum(std_db3_MC,2)./MC;

    if ESF==1
        del_sf1_Mean = sum(del_sf1_MC,2)./MC;
        del_sf2_Mean = sum(del_sf2_MC,2)./MC;
        del_sf3_Mean = sum(del_sf3_MC,2)./MC;
        std_dsf1_Mean = sum(std_dsf1_MC,2)./MC;
        std_dsf2_Mean = sum(std_dsf2_MC,2)./MC;
        std_dsf3_Mean = sum(std_dsf3_MC,2)./MC;
    end

    psi_SS_hat_Mean = sum(psi_SS_hat_MC,2)./MC;  
    std_psiSS_Mean = sum(std_psiSS_MC,2)./MC;
    the_SS_hat_Mean = sum(the_SS_hat_MC,2)./MC;  
    std_theSS_Mean = sum(std_theSS_MC,2)./MC;
    phi_SS_hat_Mean = sum(phi_SS_hat_MC,2)./MC;  
    std_phiSS_Mean = sum(std_phiSS_MC,2)./MC;

    Epsilon_Mean = sum(Epsilon_MC,2)./MC;
    std_epsilon_Mean = sum(std_epsilon_MC,2)./MC;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Specify Plotting Parameters %%

% Font size, line size, and line width %
font_size = 18;
yfont_size = 16;
leg_font_size = 12;
axes_font_size = 12;
line_size = 15;
line_width = 1;
SaveLocPNG = const.SaveLocation; % Location (filepath) to save figures to
NV = const.VecNum;               % Number of vector measurements used in EKF 
SaveLength = 75;                 % Number of runs to save figures at
SM = const.SigmaMultiplier;  	 % Value to determine how many sigma we should be plotting
plotMean = const.plotMean;       % Boolean to determine if we plot mean of monte carlo runs

if ESF == 1
    WSF = 'WithSF';
else
    WSF = 'NoSF';
end

if const.omega_S_ECI_0(3)== 0.001
    Speed = 'Slow';
end
if const.omega_S_ECI_0(3) == 0.01
    Speed = 'Spin';
end
if const.omega_S_ECI_0(3) == 0.03
    Speed = 'Fast';
end

if const.MCLength==75
    if const.inc < 10*pi/180
        SaveLocPNG = 'F:\PaperFigs\Figures\LowInc\svgImages';
        SaveLocFIG = 'F:\PaperFigs\Figures\LowInc\MATLABfigFiles';
        SaveLocDAT = 'F:\PaperFigs\Figures\LowInc\MATLABData';
    elseif const.inc > 10*pi/180 && const.inc < 60*pi/180
        SaveLocPNG = 'F:\PaperFigs\Figures\MidInc\svgImages';
        SaveLocFIG = 'F:\PaperFigs\Figures\MidInc\MATLABfigFiles';
        SaveLocDAT = 'F:\PaperFigs\Figures\MidInc\MATLABData';
    elseif const.inc > 60*pi/180
        SaveLocPNG = 'F:\PaperFigs\Figures\HighInc\svgImages';
        SaveLocFIG = 'F:\PaperFigs\Figures\HighInc\MATLABfigFiles';
        SaveLocDAT = 'F:\PaperFigs\Figures\HighInc\MATLABData';
    end
    if const.UseSGPB ==1
        Speed = 'SGPB';
        SaveLocPNG = 'F:\PaperFigs\Figures\SGPB\pngImages';
        SaveLocFIG = 'F:\PaperFigs\Figures\SGPB\MATLABfigFiles';
        SaveLocDAT = 'F:\PaperFigs\Figures\SGPB\MATLABData';
    end
end

% Create time vectors in hours for plotting %
if const.plotTime==1          % Boolean to plot in hours
    t_mag_plot = t_mag./3600;
    t_sun_plot = t_sun/3600;
    t_plot = t./3600;
    timeval='hrs';
else                          % Boolean to plot in seconds
    t_mag_plot = t_mag;
    t_sun_plot = t_sun;
    t_plot = t;
    timeval='sec';
end

% Specify plotting limits for SGPB Data, Sim Data 
if const.UseSGPB == 1
    XLimit = [0 9];
    XLimitdz = [0 1000];
else
    XLimit = [0 t_plot(end)];
    XLimitdz = [0 mu_counter];
end

plotAll = 1;
if const.MCLength == 75
    if plotAll == 1 
        MCPlot = 1:const.MCLength;
    else
        MCPlot = [5 10 15 20 25 30 35 40 45 50 55 60 65 70 75];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save Data %%
if const.MCLength==75
SaveDATName = (['75MC',Speed,WSF]); % Data Save Name
save([SaveLocDAT,'\',SaveDATName], 'Q_sv_MC','L_sv_MC','Q_sv_true_MC','L_sv_true_MC',...
    'Obs_G_sv9_MC','Obs_G_sv9_true_MC','Obs_G_sv8_MC','Obs_G_sv8_true_MC',...
    'Obs_G_sv7_MC','Obs_G_sv7_true_MC','Obs_G_sv6_MC','Obs_G_sv6_true_MC',...
    'Obs_G_sv5_MC','Obs_G_sv5_true_MC','Obs_G_sv4_MC','Obs_G_sv4_true_MC',...
    'Obs_G_sv3_MC','Obs_G_sv3_true_MC','Obs_G_sv2_MC','Obs_G_sv2_true_MC',...
    'Obs_G_sv1_MC','Obs_G_sv1_true_MC','Obs_G_cond_MC','Obs_G_cond_true_MC',...
    'Obs_G_rank_true_MC',...
    'std_epsilon_MC','Epsilon_MC','std_phiSS_MC','phi_SS_hat_MC','std_theSS_MC',...
    'the_SS_hat_MC','std_psiSS_MC','psi_SS_hat_MC','psihat_S_ECI_MC',...
    'thehat_S_ECI_MC','phihat_S_ECI_MC','b3_MC','b2_MC','b1_MC','qhat3_S_ECI_MC',...
    'qhat2_S_ECI_MC','qhat1_S_ECI_MC','qhat0_S_ECI_MC','sf3_MC','sf2_MC',...
    'sf1_MC','std_dsf3_MC','std_dsf2_MC','std_dsf1_MC','del_sf3_MC',...
    'del_sf2_MC','del_sf1_MC','std_db3_MC','std_db2_MC','std_db1_MC',...
    'del_b3_MC','del_b2_MC','del_b1_MC','std_qSS3_MC','std_qSS2_MC',...
    'std_qSS1_MC','q_SS_hat3_MC','q_SS_hat2_MC','q_SS_hat1_MC')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Relevant Parameters %%
%%%% Monte Carlo Plots %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Euler Angle Estimation Error (deg) %
PM1 = 2.5;
PM2 = 2.5;
PM3 = 2.5;
% fig=gcf;
% figure(fig.Number+1)
dEulAx(1) = subplot(311);
if plotMean ==0
    for i=MCPlot
        P2 = plot(t_plot,psi_SS_hat_MC(:,i),'r');grid on;hold on
        P1 = plot(t_plot,SM*std_psiSS_MC(:,i),'b--');
        plot(t_plot,-SM*std_psiSS_MC(:,i),'b--','HandleVisibility','off')
    end
else
    P2 = plot(t_plot,psi_SS_hat_Mean,'r');grid on;hold on
    P1 = plot(t_plot,SM*std_psiSS_Mean,'b--');
    plot(t_plot,-SM*std_psiSS_Mean,'b--','HandleVisibility','off')
end
set(gca,'FontSize',axes_font_size)
legend([P1 P2],{[num2str(SM),'\sigma'],'\delta\Psi_3_2_1 (deg)'},'FontSize',leg_font_size)
ylabel('$$\delta$$$$\psi$$','Interpreter','Latex','FontSize',yfont_size);
ylim([-2*PM1 2*PM1])
dEulAx(2) = subplot(312);
if plotMean ==0
    for i=MCPlot
        plot(t_plot,the_SS_hat_MC(:,i),'r');grid on;hold on;
        plot(t_plot,SM*std_theSS_MC(:,i),'b--');
        plot(t_plot,-SM*std_theSS_MC(:,i),'b--','HandleVisibility','off')
    end
else
    plot(t_plot,the_SS_hat_Mean,'r');grid on;hold on;
    plot(t_plot,SM*std_theSS_Mean,'b--');
    plot(t_plot,-SM*std_theSS_Mean,'b--','HandleVisibility','off') 
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$\theta$$','Interpreter','Latex','FontSize',yfont_size);
ylim([-2*PM2 2*PM2])
dEulAx(3) = subplot(313);
if plotMean ==0
    for i=MCPlot
        plot(t_plot,phi_SS_hat_MC(:,i),'r');grid on;hold on;
        plot(t_plot,SM*std_phiSS_MC(:,i),'b--');
        plot(t_plot,-SM*std_phiSS_MC(:,i),'b--','HandleVisibility','off')
    end
else
    plot(t_plot,phi_SS_hat_Mean,'r');grid on;hold on;
    plot(t_plot,SM*std_phiSS_Mean,'b--');
    plot(t_plot,-SM*std_phiSS_Mean,'b--','HandleVisibility','off')
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$\phi$$','Interpreter','Latex','FontSize',yfont_size);
ylim([-2*PM3 2*PM3])
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(dEulAx,'xy')
xlim(XLimit)
fig=gcf;
if const.MCLength==SaveLength
saveas(fig,[SaveLocPNG,'\EulErrMC',num2str(MC),'runs',WSF,Speed,'.svg'])
end

% Quaternion Estimation Error %
PM1 = 0.025;
PM2 = 0.025;
PM3 = 0.025;
fig=gcf;
figure(fig.Number+1)
dQAx(1) = subplot(311);
if plotMean ==0
    for i=MCPlot
        P2 = plot(t_plot,q_SS_hat1_MC(:,i),'r');grid on;hold on;
        P1 = plot(t_plot,SM*std_qSS1_MC(:,i),'b--');
        plot(t_plot,-SM*std_qSS1_MC(:,i),'b--','HandleVisibility','off');
    end
else
    P2 = plot(t_plot,q_SS_hat1_Mean,'r');grid on;hold on;
    P1 = plot(t_plot,SM*std_qSS1_Mean,'b--');
    plot(t_plot,-SM*std_qSS1_Mean,'b--','HandleVisibility','off'); 
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$q_{S\hat{S},1}$$','Interpreter','Latex','FontSize',font_size);
ylim([-2*PM1 2*PM1])
xlim(XLimit)
legend([P1 P2],{[num2str(SM),'\sigma'],'\deltaq'},'FontSize',leg_font_size)
dQAx(2) = subplot(312);
if plotMean ==0
    for i=MCPlot
        P2 = plot(t_plot,q_SS_hat2_MC(:,i),'r');grid on;hold on;
        P1 = plot(t_plot,SM*std_qSS2_MC(:,i),'b--');
        plot(t_plot,-SM*std_qSS2_MC(:,i),'b--','HandleVisibility','off');
    end
else
    P2 = plot(t_plot,q_SS_hat2_Mean,'r');grid on;hold on;
    P1 = plot(t_plot,SM*std_qSS2_Mean,'b--');
    plot(t_plot,-SM*std_qSS2_Mean,'b--','HandleVisibility','off'); 
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$q_{S\hat{S},2}$$','Interpreter','Latex','FontSize',font_size);
ylim([-2*PM2 2*PM2])
dQAx(3) = subplot(313);
if plotMean ==0
    for i=MCPlot
        P2 = plot(t_plot,q_SS_hat3_MC(:,i),'r');grid on;hold on;
        P1 = plot(t_plot,SM*std_qSS3_MC(:,i),'b--');
        plot(t_plot,-SM*std_qSS3_MC(:,i),'b--','HandleVisibility','off');
    end
else
    P2 = plot(t_plot,q_SS_hat3_Mean,'r');grid on;hold on;
    P1 = plot(t_plot,SM*std_qSS3_Mean,'b--');
    plot(t_plot,-SM*std_qSS3_Mean,'b--','HandleVisibility','off'); 
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$q_{S\hat{S},3}$$','Interpreter','Latex','FontSize',font_size);
ylim([-2*PM3 2*PM3])
xlim(XLimit)
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(dQAx,'xy')
fig=gcf;
if const.MCLength==SaveLength
saveas(fig,[SaveLocPNG,'\QuatErrMC',num2str(MC),'runs',WSF,Speed,'.svg'])
end

% Estimated Gyro bias Errors (deg/s) %
PM1 = 0.005;
PM2 = 0.005;
PM3 = 0.005;
fig=gcf;
figure(fig.Number+1)
dQAx(1) = subplot(311);
if plotMean ==0
    for i=MCPlot
        P2 = plot(t_plot,del_b1_MC(:,i)*r2d,'r-');grid on;hold on;
        P1 = plot(t_plot,SM*std_db1_MC(:,i)*r2d,'b--');
        plot(t_plot,-SM*std_db1_MC(:,i)*r2d,'b--','HandleVisibility','off');
    end
else
    P2 = plot(t_plot,del_b1_Mean*r2d,'r-');grid on;hold on;
    P1 = plot(t_plot,SM*std_db1_Mean*r2d,'b--');
    plot(t_plot,-SM*std_db1_Mean*r2d,'b--','HandleVisibility','off');
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$b_{g,1}$$','Interpreter','Latex','FontSize',font_size);
ylim([-2*PM1 2*PM1])
xlim(XLimit)
legend([P1 P2],{[num2str(SM),'\sigma'],'\deltab_g (deg/s)'},'FontSize',leg_font_size)
dQAx(2) = subplot(312);
if plotMean ==0
    for i=MCPlot
        P2 = plot(t_plot,del_b2_MC(:,i)*r2d,'r-');grid on;hold on;
        P1 = plot(t_plot,SM*std_db2_MC(:,i)*r2d,'b--');
        plot(t_plot,-SM*std_db2_MC(:,i)*r2d,'b--','HandleVisibility','off');
    end
else
    P2 = plot(t_plot,del_b2_Mean*r2d,'r-');grid on;hold on;
    P1 = plot(t_plot,SM*std_db2_Mean*r2d,'b--');
    plot(t_plot,-SM*std_db2_Mean*r2d,'b--','HandleVisibility','off');
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$b_{g,2}$$','Interpreter','Latex','FontSize',font_size);
ylim([-2*PM2 2*PM2])
dQAx(3) = subplot(313);
if plotMean ==0
    for i=MCPlot
        P2 = plot(t_plot,del_b3_MC(:,i)*r2d,'r-');grid on;hold on;
        P1 = plot(t_plot,SM*std_db3_MC(:,i)*r2d,'b--');
        plot(t_plot,-SM*std_db3_MC(:,i)*r2d,'b--','HandleVisibility','off');
    end
else
    P2 = plot(t_plot,del_b3_Mean*r2d,'r-');grid on;hold on;
    P1 = plot(t_plot,SM*std_db3_Mean*r2d,'b--');
    plot(t_plot,-SM*std_db3_Mean*r2d,'b--','HandleVisibility','off');
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$b_{g,3}$$','Interpreter','Latex','FontSize',font_size);
ylim([-2*PM3 2*PM3])
xlim(XLimit)
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(dQAx,'xy')
fig=gcf;
if const.MCLength==SaveLength
saveas(fig,[SaveLocPNG,'\BiasErrMC',num2str(MC),'runs',WSF,Speed,'.svg'])
end

% Estimated SF  Errors %
PM1 = 0.010;
PM2 = 0.010;
PM3 = 0.010;
fig=gcf;
figure(fig.Number+1)
dQAx(1) = subplot(311);
if plotMean ==0
    for i=MCPlot
        P2 = plot(t_plot,del_sf1_MC(:,i),'r-');grid on;hold on;
        P1 = plot(t_plot,SM*std_dsf1_MC(:,i),'b--');
        plot(t_plot,-SM*std_dsf1_MC(:,i),'b--','HandleVisibility','off');
    end
else
    P2 = plot(t_plot,del_sf1_Mean,'r-');grid on;hold on;
    P1 = plot(t_plot,SM*std_dsf1_Mean,'b--');
    plot(t_plot,-SM*std_dsf1_Mean,'b--','HandleVisibility','off');  
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$sf_{g,1}$$','Interpreter','Latex','FontSize',font_size);
if strcmp(Speed,'Fast')==1
    ylim([-2*PM1 2*PM1])
end
xlim(XLimit)
legend([P1 P2],{[num2str(SM),'\sigma'],'\deltasf_g'},'FontSize',leg_font_size)
dQAx(2) = subplot(312);
if plotMean ==0
    for i=MCPlot
        P2 = plot(t_plot,del_sf2_MC(:,i),'r-');grid on;hold on;
        P1 = plot(t_plot,SM*std_dsf2_MC(:,i),'b--');
        plot(t_plot,-SM*std_dsf2_MC(:,i),'b--','HandleVisibility','off');
    end
else
    P2 = plot(t_plot,del_sf2_Mean,'r-');grid on;hold on;
    P1 = plot(t_plot,SM*std_dsf2_Mean,'b--');
    plot(t_plot,-SM*std_dsf2_Mean,'b--','HandleVisibility','off');  
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$sf_{g,2}$$','Interpreter','Latex','FontSize',font_size);
if strcmp(Speed,'Fast')==1
    ylim([-2*PM1 2*PM1])
end
xlim(XLimit)
dQAx(3) = subplot(313);
if plotMean ==0
    for i=MCPlot
        P2 = plot(t_plot,del_sf3_MC(:,i),'r-');grid on;hold on;
        P1 = plot(t_plot,SM*std_dsf3_MC(:,i),'b--');
        plot(t_plot,-SM*std_dsf3_MC(:,i),'b--','HandleVisibility','off');
    end
else
    P2 = plot(t_plot,del_sf3_Mean,'r-');grid on;hold on;
    P1 = plot(t_plot,SM*std_dsf3_Mean,'b--');
    plot(t_plot,-SM*std_dsf3_Mean,'b--','HandleVisibility','off');  
end
set(gca,'FontSize',axes_font_size)
ylabel('$$\delta$$$$sf_{g,3}$$','Interpreter','Latex','FontSize',font_size);
if strcmp(Speed,'Fast')==1
    ylim([-2*PM1 2*PM1])
end
xlim(XLimit)
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(dQAx,'xy')
xlim(XLimit)
fig=gcf;
if const.MCLength==SaveLength
saveas(fig,[SaveLocPNG,'\SFErrMC',num2str(MC),'runs',WSF,Speed,'.svg'])
end

% Angular Error between Estimate and True Spacecraft 3 axis
fig=gcf;
figure(fig.Number+1)
if plotMean ==0
    STDSum=0;
    for i=MCPlot
        P2 = plot(t_plot,Epsilon_MC(:,i),'r');hold on;
        STDSum = STDSum+std_epsilon_MC(:,i);
    end
    P1 = plot(t_plot,STDSum./MC*SM,'b--');
else
    P2 = plot(t_plot,Epsilon_Mean,'r');hold on;
    P1 = plot(t_plot,std_epsilon_Mean*SM,'b--');
end
set(gca,'FontSize',axes_font_size)
grid on; ylabel('$$\epsilon$$','Interpreter','Latex','FontSize',font_size);
legend([P1 P2],{[num2str(SM),'\sigma'],'\epsilon (deg)'},'FontSize',leg_font_size)
% yline(const.KnowledgeError,...
%     'k--',['Knowledge Error Requirement = ',...
%     num2str(const.KnowledgeError),' deg'],'FontSize',leg_font_size,'LineWidth',1,'HandleVisibility','off');
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
xlim(XLimit)
ylim([0 10])
fig=gcf;
if const.MCLength==SaveLength
saveas(fig,[SaveLocPNG,'\EpsilonMC',num2str(MC),'runs',WSF,Speed,'.svg'])
end

if const.RunSOTest ==1
    % Deterministic Observability % 
    len_plot = 10;
    fig=gcf;
    figure(fig.Number+1)
    if plotMean ==0
        for i=MCPlot
            scatter(t_mag,Obs_G_rank_MC(:,i),'r*');grid on; hold on;
            scatter(t_mag,Obs_G_rank_true_MC(:,i),'k.')
        end
    else
        scatter(t_mag,Obs_G_rank_Mean,'r*');grid on; hold on;
        scatter(t_mag,Obs_G_rank_true_Mean,'k.')
    end
    ylabel('Rank \it $$O_{k}^{T}O_{k}$$','Interpreter','Latex','FontSize',font_size);
    xlabel(['Time (sec)'],'Interpreter','Latex','FontSize',font_size);
    legend('Noisy System','True System')
    title('Deterministic Observability','Interpreter','Latex','FontSize',font_size)
    ylim([0 FullRank+2])
    xlim([0 len_plot*1/const.MagSamp])
    fig=gcf;
    if const.MCLength==SaveLength
    saveas(fig,[SaveLocPNG,'\DetermObsv',num2str(MC),'runs',WSF,Speed,'.svg'])
    end

    % Stochastic Observability %
    fig=gcf;
    figure(fig.Number+1)
    subplot(211);
    if plotMean ==0
        for i=MCPlot
            plot(t_mag,L_sv_MC(:,i),'r-');grid on; hold on;
            plot(t_mag,L_sv_true_MC(:,i),'k--');
        end
    else
        plot(t_mag,L_sv_Mean,'r-');grid on; hold on;
        plot(t_mag,L_sv_true_Mean,'k--');     
    end
    ylabel('$$\sigma_{max}$$($$\Lambda_{k})$$','Interpreter','Latex','FontSize',font_size);
    xlim([0 len_plot*1/const.MagSamp])
    xlabel(['Time (sec)'],'Interpreter','Latex','FontSize',font_size);
    legend('Noisy System','True System')
    title('Stochastic Observability','Interpreter','Latex','FontSize',font_size)
    subplot(212);
    if plotMean ==0
        for i=MCPlot
            plot(t_mag,Q_sv_MC(:,i),'r-');grid on; hold on;
            plot(t_mag,Q_sv_true_MC(:,i),'k--');
        end
    else
        plot(t_mag,Q_sv_Mean,'r-');grid on; hold on;
        plot(t_mag,Q_sv_true_Mean,'k--');     
    end
    ylabel('$$\sigma_{max}$$($$\overline {Q}_{k})$$','Interpreter','Latex','FontSize',font_size);
    ylim([0 5e-3])
    xlabel(['Time (hrs)'],'Interpreter','Latex','FontSize',font_size);
    fig=gcf;
    if const.MCLength==SaveLength
    saveas(fig,[SaveLocPNG,'\StochObsv',num2str(MC),'runs',WSF,Speed,'.svg'])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
