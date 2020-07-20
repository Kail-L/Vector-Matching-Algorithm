%% Vector Matching Algorithm Extended Kalman Filter (EKF) for CubeSat %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to solve attitude estimation problem for sattelite truth data
% simulated from SatAttEst.m.
%
% Kail Laughlin
% 7/9/20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set quadprog options
optionsQP =  optimset('Algorithm','active-set','Display','off'); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Constants & Set Up Data %%

% Data from post_processing.m
%       IMU = [t, Gyro_Good, Gyro_Corrupt]
%       TrueAng = [t psi psidot theta thetadot phi phidot]
%       Mag = [B_ECI B_S]
%       Sun = [S_ECI S_S S_S_noisy]

% Constants %
constants_struct
I_b = const.I_b;                    % CubeSat Moment of Inertia
DoReboot = const.PerformDailyReboot;% Boolean to say if we do a daily reboot or not
DR = const.DailyRebootTime;         % Daily Reboot Time, hr
tau = const.G_tau;                  % gyro correlation time for gauss markov, s
IRb = const.G_inRunBiasSTD*d2r/3600;% In-Run bias, rad/s
ARW = (const.G_ARW/60*d2r);         % rad/sqrt(sec)

% Determine which noise to populate Q_w with %
if const.G_UseARW == 0
    OutNoise = const.G_measNoise*d2r;   % rad/s
else
    OutNoise = ARW;
end

UseMeas = const.UseMeas;            % Boolean to say if we should perform measurements
Magvar = const.M_measNoise;         % Magnetometer variance, mgauss
Sunvar = const.S_measNoise*d2r;     % Sun sensor variance, rad
ASF = const.G_AddScaleFactor;       % Boolean used to add scale factor to gyro measurements
ESF = const.EstimateSF;             % Boolean used to estimate scale factor
SFnoise = const.G_SFNoise;          % Scale factor noise, ppm
SFLB=const.G_ScaleFactorLB;         % Lower bound scale factor matrix estimate for quadprog
SFUB=const.G_ScaleFactorUB;         % Upper bound on scale factor matrix estimate for quadprog
IFlag = const.IQ;                   % Flag to determine quaternion initialization
SR = const.MagSamp;                 % Measurement update rate
CSg = const.CSg;                    % Misalignment DCM between gyro and body 
RebootLength = const.RebootLength;  % Time that System is undergoing reboot
RebootOccurred = 0;                 % Flag to determine if reboot occured
alpha = const.alpha;                % Tuning Parameter for G matrix
beta = const.beta;                  % Tuning Parameter for G matrix
gamma = const.gamma;                % Tuning parameter for G matrix

% NOTE: The tuning parameters above are used to ensure our estimate
% converges in an acceptable amount of time. However, the larger the tuning
% parameters, the larger the covariance bound will be on the estimate.
% Typically, there is a trade off between convergence speed and accuracy.

% Data %
t = IMU(:,1);                       % time vector, s
dt=mean(diff(t));                   % Time step of gyro timing data from IMU_SV_Mag.mat
drl=length(t);                      % Looping variable for integration

uECI=Mag(:,1:3)*T2mg;                % Mag field vector in inertial frame, mGauss
Mag_Clean_S=Mag(:,4:6)*T2mg;         % Clean Mag field vector in body frame, mGauss 

sECI=Sun(:,1:3);                    % Unit sun vector measurement in inertial frame
Sun_Clean_S=Sun(:,4:6);             % Unit sun vector measurement in body frame
sS=Sun(:,7:9);                      % Unit sun vector measurment in body frame

Om_SE_good=IMU(:,2:4);              % Truth angular rate values 
psi_S_ECI=TrueAng(:,2)*r2d;         % Truth value of psi, deg
the_S_ECI=TrueAng(:,4)*r2d;         % Truth value of theta, deg
phi_S_ECI=TrueAng(:,6)*r2d;         % Truth value of phi, deg
q_S_ECI=q_S_ECI_true;               % True quaternion (note, this is unecessary but placed here for clarity)

% Initialize Matrices %
qhat_S_ECI = zeros(drl,4);          % Estimate quaternion
q_SS_hat = zeros(drl,4);            % Quaternion errors from body frame to estimate body frame
psihat_S_ECI=zeros(drl,1);          % Estimate value of psi, deg
thehat_S_ECI=zeros(drl,1);          % Estimate value of theta, deg
phihat_S_ECI=zeros(drl,1);          % Estimate value of phi, deg
psi_SS_hat = zeros(drl,1);          % Error in psi, deg
the_SS_hat = zeros(drl,1);          % Error in theta, deg
phi_SS_hat = zeros(drl,1);          % Error in phi, deg


del_q_SS = zeros(drl,3);            % Quaternion vector estimation errors
del_b = zeros(drl,3);               % Bias estimation errors
del_sf = zeros(drl,3);              % Scale factor estimation errors

b = zeros(drl,3);                   % Storage vector of biases
b_out = [0 0 0]'*d2r;               % Initial bias estimate on gyros
b(1,:) = b_out';                    % Initial bias estimate
sf = zeros(drl,3);                  % Storage vector of scale factors
sf_out = [0 0 0]';                  % Initial scale factor estimate
sf(1,:) = sf_out';                  % Initial scale factor estimate

dz_plot = zeros(length(t_mag),6);   % Innovations error
Epsilon = zeros(drl,1);             % Angular error between estimate and true body 3 axes
MagAng = zeros(length(t_mag),1);    % Angle between initial magnetic vector and current magnetic vector
delMagAng = zeros(length(t_mag),1); % Angle between current and prior EMF vectors
uS_rec = zeros(length(t_mag),3);    % Magnetic field in body frame measured                              
uS_comp = zeros(length(t_mag),3);   % Magnetic field in body frame computed with true attitude
sS_rec = zeros(length(t_sun),3);    % Sun vector in body frame measured
sS_comp = zeros(length(t_sun),3);   % Sun vector in body frame computed with true attitude

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% State Covariance and Storage %%
if ESF ==1
    P_stor = zeros(drl,9);          % Covariance Storage vector
else
    P_stor = zeros(drl,6);
end
P_eul = zeros(4,4);                 % Covariance matrix used to convert quat covariance to eul covariance

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
%% Observability Values for VMA_SO_Test.m %%
if ESF ==1
    Fd_TM = eye(9);                 % State transition matrix of F
else
    Fd_TM = eye(6);
end

Obs = [];                           % Observability Grammian
Obs_rank = zeros(length(t_mag),1);  % Storage vector for rank of observability matrix
Obs_G_rank = zeros(length(t_mag),1);% Observability grammian rank check

% Stochastic Observability Vals %
H_Test = cell(length(t_mag),1);       % Storage cell of measurement matrices
H_Test_true = cell(length(t_mag),1);  % Storage cell of measurement matrices 
H_row = 0;                          % Initial value of number of rows of H
H_col = 0;                          % Initial value of number of columns of H
PHI_Test = cell(drl,1);               % Storage cell for current state transition matrix
PHI_Test_true = cell(drl,1);          % Storage cell for current state transition matrix
PHI_Test{1} = Fd_TM;                  % Initial state transistion matrix
PHI_Test_true{1} = Fd_TM;             % Initial state transistion matrix
PhiSize = size(Fd_TM);              % Size of state transition matrix
PHI_row = PhiSize(1);               % Rows of state transition matrix
PHI_col = PhiSize(2);               % Columns of state transition matrix
Rz_Test = cell(drl,1);                % Storage cell for measurement covariance matrix
Rz_Test_true = cell(drl,1);           % Storage cell for measurement covariance matrix
Qd_Test = cell(drl,1);                % Storage cell for process noise covariance matrix
Qd_Test_true = cell(drl,1);           % Storage cell for process noise covariance matrix
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
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corrupt Gyro Parameters %
[Om_SE_noisy,G_SF,B_true] = emulateGyro(t,Om_SE_good,const);

% Difference in True and Noisy Angular Velocity
Om_SE_difference = (Om_SE_noisy-Om_SE_good);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Corrupt Mag Parameters %
if const.UseSGPB == 1
    uS = B_S*T2mg;      % Convert to mGauss
else
    uS = emulateMag(t_mag,Mag_Clean_S,const);
%     uS = CalMag;
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
    std_dsf1(1) = sqrt(P0(7,7));
    std_dsf2(1) = sqrt(P0(8,8));
    std_dsf3(1) = sqrt(P0(9,9));
else % If we are NOT estimating scale factor
    P_stor(1,:) = [P0(1,1) P0(2,2) P0(3,3) P0(4,4) P0(5,5) P0(6,6)];                
    std_qSS1(1) = sqrt(P0(1,1));
    std_qSS2(1) = sqrt(P0(2,2));
    std_qSS3(1) = sqrt(P0(3,3));
    std_db1(1) = sqrt(P0(4,4));
    std_db2(1) = sqrt(P0(5,5));
    std_db3(1) = sqrt(P0(6,6));
end

state_cov_his{1}=P0;           % Store for SO test

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sensor Noise Covariance %
Rh = ((Magvar)^2)*eye(3);         % Magnetometer covariance matrix, mGauss
Rs = ((Sunvar)^2)*eye(3);         % Sun Sensor covariance matrix
if const.UseSun==1
    Rz(1:3,1:3) = Rh;             % State Measurement covariance matrix
    Rz(4:6,4:6) = Rs;
else
    Rz = Rh;                      % State Measurement covariance matrix
end

Rz_Test{1} = Rz;                    % Store for SO test
Rz_Test_true{1} = Rz;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process Noise Covariance %
if ESF ==1
    Qw = [(OutNoise)^2*eye(3) zeros(3,3) zeros(3,3);    % Gyro Output noise
          zeros(3,3) 2*(IRb)^2/tau*eye(3) zeros(3,3);   % Gyro in-run bias repeat.
          zeros(3,3) zeros(3,3) (SFnoise/1e6)^2*eye(3)];% Scale factor noise 
else
    Qw = [(OutNoise)^2*eye(3) zeros(3,3);   % Gyro Output noise
          zeros(3,3) 2*(IRb)^2/tau*eye(3)]; % Gyro in-run bias repeat.
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Quaternion %%
q_S_ECI_0 = q_S_ECI(1,:);           % Initial truth quaternion

if IFlag == 0                       % IFlag == 0: Use truth quat as initial estimate
    qhat_S_ECI_0 = q_S_ECI_0';  
elseif IFlag == 1                   % IFlag == 1: Use max angular error allowed, deg
    MaxAngErr=const.ME;             % Maximum angular error allowed, deg
    ErrEul = [MaxAngErr*(2*rand-1);
              MaxAngErr*(2*rand-1);
              MaxAngErr*(2*rand-1)];
    ErrEul = ErrEul*d2r;
    C_SECIErr=eul2dcm(ErrEul);
    Chat_S_ECI_0 = C_SECIErr*Cba_0;
    qhat_S_ECI_0 = dcm2quat(Chat_S_ECI_0)';
elseif IFlag == 2                   % IFlag == 2: Use completely random quat          
    q0=rand(1,1);                   % Random quaternion
    q1=2*rand(1,1) - 1;
    q2=2*rand(1,1) - 1;
    q3=2*rand(1,1) - 1;
    qhat_S_ECI_0 = [q0 q1 q2 q3]';
    qhat_S_ECI_0 = qhat_S_ECI_0/quatnorm(qhat_S_ECI_0);
elseif IFlag == 3                   % IFlag == 3: Use set of error Euler angles
    ErrEul=const.ErrEul;            % Euler angle error
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
psihat_S_ECI_0 = EE_hat_0(1);   % yaw
thehat_S_ECI_0 = EE_hat_0(2);   % pitch
phihat_S_ECI_0 = EE_hat_0(3);   % roll

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Determine Initial Angular Error, Epsilon %%
Chat_S_ECI_0 = quat2dcm(qhat_S_ECI_0');
C_SS_hat_0 = quat2dcm(q_S_ECI_0)*inv(Chat_S_ECI_0);
vS = [0 0 1]';                  % 3-axis of S frame
vS_hat = C_SS_hat_0*vS;         % 3-axis of S_hat frame
epsilon = acos(C_SS_hat_0(3,3))*180/pi;

disp(['EKF Initial Yaw Estimate:   ',num2str(psihat_S_ECI_0),' deg'])
disp(['EKF Initial Pitch Estimate: ',num2str(thehat_S_ECI_0),' deg'])
disp(['EKF Initial Roll Estimate:  ',num2str(phihat_S_ECI_0),' deg'])
fprintf('\n')
disp(['Initial Angular Error Between S3 and S3_hat: ',num2str(epsilon),' degrees'])
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize quaternion error, body frame %%
q_SS_hat_0 = quatmult(q_S_ECI_0',quatInverse(qhat_S_ECI_0'));
q_SS_hat_0 = q_SS_hat_0/norm(q_SS_hat_0);
q_SS_hat(1,:) = q_SS_hat_0';
del_q_SS(1,:) = q_SS_hat(1,2:4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Yaw, Pitch, Roll, degrees %
eulhat_ba_0 = dcm2eul(Chat_S_ECI_0);
psihat_S_ECI(1) = eulhat_ba_0(1)*r2d;
thehat_S_ECI(1) = eulhat_ba_0(2)*r2d;
phihat_S_ECI(1) = eulhat_ba_0(3)*r2d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run EKF %%
mu_counter=2;           % Measurement update counter 
disp('Running EKF...')
fprintf('\n')
for k = 2:drl
%%% Check If Daily Reboot Occured %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if DoReboot==1 && t(k)>=(DR*3600) && t(k)<=((DR*3600)+RebootLength)
        RebootOccurred=1;
        P0=Pinit;
        std_psiSS(k) = 100;     % NOTE: no way to calculate this val without quaternion error, so assume very large (deg)
        std_theSS(k) = 100;     % NOTE: no way to calculate this val without quaternion error, so assume very large (deg)
        std_phiSS(k) = 100;     % NOTE: no way to calculate this val without quaternion error, so assume very large (deg)
        std_qSS1(k) = sqrt(P0(1,1));
        std_qSS2(k) = sqrt(P0(2,2));
        std_qSS3(k) = sqrt(P0(3,3));
        std_db1(k) = sqrt(P0(4,4));
        std_db2(k) = sqrt(P0(5,5));
        std_db3(k) = sqrt(P0(6,6));
        if ASF ==1
            std_dsf1(k) = sqrt(P0(7,7));
            std_dsf2(k) = sqrt(P0(8,8));
            std_dsf3(k) = sqrt(P0(9,9));
        end
        if mu_counter<=length(t_mag)&&t_mag(mu_counter)==t(k)
            mu_counter=mu_counter+1;
        end
        if t(k)==(DR*3600)
            q_prior = qhat_S_ECI(k-1,:);
            w_prior = Om_SE_noisy(k-1,:)'-b_out;
        end
        if t(k)==((DR*3600)+RebootLength)
            [qhat_ba_rebooted,what_ba_rebooted,RTime]=rebootAtt(q_prior,w_prior,1000,const);
            qhat_S_ECI(k,:)=qhat_ba_rebooted;
        end
    else
%%%%%%%% Time Update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ESF ==1  % If we are estimating scale factors
            %Estimated Trajectory%
            iSFmat = inv(diag([1+sf_out(1) 1+sf_out(2) 1+sf_out(3)]));
            Om_hat_SE = iSFmat*CSg*(Om_SE_noisy(k,:)'-b_out);% Estimate angular velocity from S frame to ECI frame
            
            F = [-sk(Om_hat_SE) -1/2*(iSFmat)*CSg -1/2*diag(CSg*(Om_SE_noisy(k,:)'-b_out)); 
                zeros(3) (-1/tau)*eye(3) zeros(3);
                zeros(3) zeros(3) zeros(3)];

            G = [-1/2*iSFmat*alpha*CSg zeros(3) zeros(3);         
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
        del_b(k,:) = del_b(k-1,:);
        del_sf(k,:) = del_sf(k-1,:);
        del_q_SS(k,:) = del_q_SS(k-1,:);
             
        % Covariance Update 
        Phi = expm(F*dt);           % Discreet time F matrix (Phi matrix)
        PhiTrue = expm(Ftrue*dt);   % Discreet time Ftrue matrix
        Cd = disrw(F,G,dt,Qw);      % Discreet time noise mapping matrix (Gamma matrix)
        CdTrue = disrw(Ftrue,Gtrue,dt,Qw);% Discreet time noise mapping matrix
        P = Phi*P0*Phi' + Cd;       % Updated covariance
        PhiSize = size(Phi);
        CdSize = size(Cd);        
                
        % Estimate Quaternion Update 
        qhat_S_ECI(k,:) = quatIntegrate((Om_hat_SE)',qhat_S_ECI(k-1,:)',dt);
        qhat_S_ECI(k,:) = qhat_S_ECI(k,:)/quatnorm(qhat_S_ECI(k,:));
        if (qhat_S_ECI(k,1) < 0)         % Enforce Positive Scalar 
            qhat_S_ECI(k,:) = -qhat_S_ECI(k,:);
        end
 
        % Enforce Positive Scalar on True Quaternion 
        if q_S_ECI(k,1) < 0              % Enforce Positive Scalar 
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
            
            % Current Estimate DCM
            Chat_S_ECI = quat2dcm(qhat_S_ECI(k,:));
            C_S_ECI = quat2dcm(q_S_ECI(k,:));
            
            % Non-Normalized Mag Field model & measurement 
            US = uS(mu_counter,:)';     uS_rec(mu_counter,:) = US';
            UECI = uECI(mu_counter,:)';
            
            % Normalized Sun Sensor model & measurement (if used) %
            SS = sS(mu_counter,:)';     sS_rec(mu_counter,:) = SS';
            SECI = sECI(mu_counter,:)';  
            
            % Current EMF vector, body frame %
            MagVecCurrent = US/norm(US);
            
            % Truth values of mag field & sun vec in body frame %
            uS_comp(mu_counter,:) = (eul2dcm([psi_S_ECI(k) the_S_ECI(k) phi_S_ECI(k)]*d2r)*UECI)';
            sS_comp(mu_counter,:) = (eul2dcm([psi_S_ECI(k) the_S_ECI(k) phi_S_ECI(k)]*d2r)*SECI)';
            
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
            
            % Store values for stochastic observability test
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
                dz_plot(mu_counter,:) = [dmag' dsun'];
            else 
                dz = dmag;
                dz_plot(mu_counter,:) = [dmag' dsun'];
            end
            
%%%%%%%%%%%% Run EKF%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
            Rz_Test{mu_counter} = Rz;  
            Rz_Test_true{mu_counter} = Rz;
            SizeP = size(P);
            K = (P*H')/(H*P*H' + Rz); % Kalman Filter Gain 
            P = (eye(SizeP(1)) - K*H)*P;      % Covariance update 
            del_x = K*dz;                     % State Error Vector 

%%% Constrain SF estimate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
             if ESF ==1
                if const.constrainEst == 1
                % Constrain using QuadProg
                    D=[];
                    d=[];
                    lbqp = [];
                    ubqp = [];
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
            if ESF ==1
                del_q_SS(k,:) = xbar(1:3)';
                del_b(k,:) = xbar(4:6)';   
                del_sf(k,:) = xbar(7:9);
            else
                del_q_SS(k,:) = del_x(1:3)';
                del_b(k,:) = del_x(4:6)';
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
        q_SS_hat(k,:) = quatmult(q_S_ECI(k,:)',quatInverse(qhat_S_ECI(k,:)));
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
        std_epsilon(k) = sqrt(Peps);
    end
end 

% Display Values to Command Line %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ScaleFactor = [sf_out(1) sf_out(2) sf_out(3)];
Bias = b_out*r2d;

disp(['Final Gyro Bias Estimate, S1 axis: ',num2str(Bias(1)),' deg/s'])
disp(['Final Gyro Bias Estimate, S2 axis: ',num2str(Bias(2)),' deg/s'])
disp(['Final Gyro Bias Estimate, S3 axis: ',num2str(Bias(3)),' deg/s'])
disp(['Final Gyro SF Estimate, S1 axis: ',num2str(ScaleFactor(1))])
disp(['Final Gyro SF Estimate, S2 axis: ',num2str(ScaleFactor(2))])
disp(['Final Gyro SF Estimate, S3 axis: ',num2str(ScaleFactor(3))])
if RebootOccurred == 1
    disp(['Attitude Reboot Computation Time: ',num2str(RTime),' s'])
end
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Stochastic Observability Test %%
if const.RunSOTest ==1
    disp('Running Observability Test...')
    fprintf('\n')
    
    % Estimate Trajectory %
    PHI_SO = PHI_Test;
    Qd_SO = Qd_Test;
    H_SO = H_Test;
    Rz_SO = Rz_Test;
    VMA_SO_Test                         % Run SO test
    Obs_G_rank = r_OBS_gram';           % Rank of observability gramian
    Obs_G_cond = c_OBS_gram';           % Condition # of observability gramian
    Obs_G_sv1 = OBS_gram_sv(1,:)';      % 1st singular value of obsv. gram.
    Obs_G_sv2 = OBS_gram_sv(2,:)';      % 2nd singular value of obsv. gram.
    Obs_G_sv3 = OBS_gram_sv(3,:)';      % 3rd singular value of obsv. gram.
    Obs_G_sv4 = OBS_gram_sv(4,:)';      % 4th singular value of obsv. gram.
    Obs_G_sv5 = OBS_gram_sv(4,:)';      % 5th singular value of obsv. gram.
    Obs_G_sv6 = OBS_gram_sv(6,:)';      % 6th singular value of obsv. gram.
    if ESF==1
        Obs_G_sv7 = OBS_gram_sv(7,:)';  % 7th singular value of obsv. gram.
        Obs_G_sv8 = OBS_gram_sv(7,:)';  % 8th singular value of obsv. gram.
        Obs_G_sv9 = OBS_gram_sv(7,:)';  % 9th singular value of obsv. gram.
    end
    L_sv = (Lambda_sv(1,:))';      % Max singular value of Lambda metric (see Chapter 4. of Kail's Thesis)
    Q_sv = (Q_bar_sv(1,:))';        % Max singular value of Q bar metric (see Chapter 4. of Kail's Thesis)

    % Truth Trajectory %
    PHI_SO = PHI_Test_true;
    Qd_SO = Qd_Test_true;
    H_SO = H_Test_true;
    Rz_SO = Rz_Test_true;
    VMA_SO_Test                              % Run SO test
    Obs_G_rank_true = r_OBS_gram';           % Rank of observability gramian
    Obs_G_cond_true = c_OBS_gram';           % Condition # of observability gramian
    Obs_G_sv1_true = OBS_gram_sv(1,:)';      % 1st singular value of obsv. gram.
    Obs_G_sv2_true = OBS_gram_sv(2,:)';      % 2nd singular value of obsv. gram.
    Obs_G_sv3_true = OBS_gram_sv(3,:)';      % 3rd singular value of obsv. gram.
    Obs_G_sv4_true = OBS_gram_sv(4,:)';      % 4th singular value of obsv. gram.
    Obs_G_sv5_true = OBS_gram_sv(4,:)';      % 5th singular value of obsv. gram.
    Obs_G_sv6_true = OBS_gram_sv(6,:)';      % 6th singular value of obsv. gram.
    if ESF==1
        Obs_G_sv7_true = OBS_gram_sv(7,:)';  % 7th singular value of obsv. gram.
        Obs_G_sv8_true = OBS_gram_sv(7,:)';  % 8th singular value of obsv. gram.
        Obs_G_sv9_true = OBS_gram_sv(7,:)';  % 9th singular value of obsv. gram.
    end
    L_sv_true = (Lambda_sv(1,:))';       % Max singular value of Lambda metric (see Chapter 4. of Kail's Thesis)
    Q_sv_true = (Q_bar_sv(1,:))';        % Max singular value of Q bar metric (see Chapter 4. of Kail's Thesis)
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
NV = const.VecNum;            % Number of vector measurements used in EKF
SM = const.SigmaMultiplier;   % Value to determine how many sigma we should be plotting

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot the Output of Relevant Parameters %%

% Clean and Noisy Omega %
subplot(321);
plot(t,r2d*Om_SE_good(:,1),'r');
grid;ylabel('$\omega^{SE}_{S,1}$','fontsize',font_size,'Interpreter','latex');
title('Error Free $\omega$ (deg/s)','fontsize',font_size,'Interpreter','latex')
subplot(322);
plot(t,r2d*Om_SE_noisy(:,1),'g');
grid;ylabel('$\omega^{SE}_{S,1}$','fontsize',font_size,'Interpreter','latex');
title('Corrupted $\omega$ (deg/s)','fontsize',font_size,'Interpreter','latex')
subplot(323);
plot(t,r2d*Om_SE_good(:,2),'r');grid;
ylabel('$\omega^{SE}_{S,2}$','fontsize',font_size,'Interpreter','latex');
subplot(324);
plot(t,r2d*Om_SE_noisy(:,2),'g');
grid;ylabel('$\omega^{SE}_{S,2}$','fontsize',font_size,'Interpreter','latex');
subplot(325);
plot(t,r2d*Om_SE_good(:,3),'r');
grid;ylabel('$\omega^{SE}_{S,3}$','fontsize',font_size,'Interpreter','latex');
xlabel('Time (sec)','fontsize',font_size,'Interpreter','latex');
subplot(326)
plot(t,r2d*Om_SE_noisy(:,3),'g');
grid;ylabel('$\omega^{SE}_{S,3}$','fontsize',font_size,'Interpreter','latex');
xlabel('Time (sec)','fontsize',font_size,'Interpreter','latex');

% Clean and Noisy EMF %
fig=gcf;
figure(fig.Number+1)
subplot(321);
plot(t_mag,Mag_Clean_S(:,1),'r');
grid;ylabel('$B^{S}_{S,1}$','fontsize',font_size,'Interpreter','latex');
title('Error Free $B$ (mG)','fontsize',font_size,'Interpreter','latex')
subplot(322);
plot(t_mag,uS(:,1),'g');
grid;
ylabel('$B^{S}_{S,1}$','fontsize',font_size,'Interpreter','latex');
title('Corrupted $B$ (mG)','fontsize',font_size,'Interpreter','latex')
subplot(323);
plot(t_mag,Mag_Clean_S(:,2),'r');grid;
ylabel('$B^{S}_{S,2}$','fontsize',font_size,'Interpreter','latex');
subplot(324);
plot(t_mag,uS(:,2),'g');grid;
ylabel('$B^{S}_{S,2}$','fontsize',font_size,'Interpreter','latex');
subplot(325);
plot(t_mag,Mag_Clean_S(:,3),'r');
grid;
ylabel('$B^{S}_{S,3}$','fontsize',font_size,'Interpreter','latex');
xlabel('Time (sec)','fontsize',font_size,'Interpreter','latex');
subplot(326)
plot(t_mag,uS(:,3),'g');
grid;ylabel('$B^{S}_{S,3}$','fontsize',font_size,'Interpreter','latex');
xlabel('Time (sec)','fontsize',font_size,'Interpreter','latex');

% Clean and Noisy Sun Vector %
if const.UseSun ==1
    fig=gcf;
    figure(fig.Number+1)
    subplot(321);
    plot(t_sun,Sun_Clean_S(:,1),'r');
    grid;ylabel('$S^{S}_{S,1}$','fontsize',font_size,'Interpreter','latex');
    title('Error Free $S$ Vector','fontsize',font_size,'Interpreter','latex')
    subplot(322);
    plot(t_sun,sS(:,1),'g');
    grid;
    ylabel('$S^{S}_{S,1}$','fontsize',font_size,'Interpreter','latex');
    title('Corrupted $S$ Vector','fontsize',font_size,'Interpreter','latex')
    subplot(323);
    plot(t_sun,Sun_Clean_S(:,2),'r');grid;
    ylabel('$S^{S}_{S,2}$','fontsize',font_size,'Interpreter','latex');
    subplot(324);
    plot(t_sun,sS(:,2),'g');grid;
    ylabel('$S^{S}_{S,2}$','fontsize',font_size,'Interpreter','latex');
    subplot(325);
    plot(t_sun,Sun_Clean_S(:,3),'r');
    grid;
    ylabel('$S^{S}_{S,3}$','fontsize',font_size,'Interpreter','latex');
    xlabel('Time (sec)','fontsize',font_size,'Interpreter','latex');
    subplot(326)
    plot(t_sun,sS(:,3),'g');
    grid;ylabel('$S^{S}_{S,3}$','fontsize',font_size,'Interpreter','latex');
    xlabel('Time (sec)','fontsize',font_size,'Interpreter','latex');
end

% 3D plot of True EMF in ECI frame %
fig=gcf;
figure(fig.Number+1)
plot3(uECI(:,1),uECI(:,2),uECI(:,3),'ro')
xlabel('$$B^{E}_{E,1}$$','Interpreter','Latex','FontSize',font_size);
ylabel('$$B^{E}_{E,2}$$','Interpreter','Latex','FontSize',font_size);
zlabel('$$B^{E}_{E,3}$$','Interpreter','Latex','FontSize',font_size);
grid on; 
title('Error Free $B$ (mG)','fontsize',font_size,'Interpreter','latex')

% 3D plot of True EMF in S frame %
fig=gcf;
figure(fig.Number+1)
plot3(Mag_Clean_S(:,1),Mag_Clean_S(:,2),Mag_Clean_S(:,3),'bo')
xlabel('$$B^{S}_{S,1}$$','Interpreter','Latex','FontSize',font_size);
ylabel('$$B^{S}_{S,2}$$','Interpreter','Latex','FontSize',font_size);
zlabel('$$B^{S}_{S,3}$$','Interpreter','Latex','FontSize',font_size);
grid;
title('Error Free $B$ (mG)','fontsize',font_size,'Interpreter','latex')

% 3D plot of Noisy EMF in S frame %
fig=gcf;
figure(fig.Number+1)
plot3(uS(:,1),uS(:,2),uS(:,3),'go')
xlabel('$$B^{S}_{S,1,meas}$$','Interpreter','Latex','FontSize',font_size);
ylabel('$$B^{S}_{S,2,meas}$$','Interpreter','Latex','FontSize',font_size);
zlabel('$$B^{S}_{S,3,meas}$$','Interpreter','Latex','FontSize',font_size);
grid;
title('Corrupted $B$ (mG)','fontsize',font_size,'Interpreter','latex')

% True & Estimated Euler Angles %
fig=gcf;
figure(fig.Number+1)
EulAx(1) = subplot(311);
plot(t_plot,psihat_S_ECI,'r--',t_plot,psi_S_ECI,'b--');grid on;
ylabel('$$\psi$$','Interpreter','Latex','FontSize',font_size);
title('Euler Angles (deg)','Interpreter','Latex','FontSize',font_size);
legend('Estimate','True','Location','Best')
EulAx(2) = subplot(312);
plot(t_plot,thehat_S_ECI,'r--',t_plot,the_S_ECI,'b--');grid on;
ylabel('$$\theta$$','Interpreter','Latex','FontSize',font_size);
EulAx(3) = subplot(313);
plot(t_plot,phihat_S_ECI,'r--',t_plot,phi_S_ECI,'b--');grid on;
ylabel('$$\phi$$','Interpreter','Latex','FontSize',font_size);
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(EulAx,'xy')
xlim(XLimit)

% Euler Angle Estimation Error %
PM1 = 4;
PM2 = 4;
PM3 = 4;
fig=gcf;
figure(fig.Number+1)
dEulAx(1) = subplot(311);
plot(t_plot,psi_SS_hat,'r-');grid on;hold on;
ylabel('$$\delta$$$$\psi$$','Interpreter','Latex','FontSize',yfont_size);
plot(t_plot,SM*std_psiSS,'b')
plot(t_plot,-SM*std_psiSS,'b')
ylim([-2*PM1 2*PM1])
legend('Error',[num2str(SM),'\sigma STD'])
title('Euler Angle Estimation Error (deg)','Interpreter','Latex','FontSize',font_size);
dEulAx(2) = subplot(312);
plot(t_plot,the_SS_hat,'r-');grid on;hold on;
ylabel('$$\delta$$$$\theta$$','Interpreter','Latex','FontSize',yfont_size);
plot(t_plot,SM*std_theSS,'b')
plot(t_plot,-SM*std_theSS,'b')
ylim([-2*PM2 2*PM2])
dEulAx(3) = subplot(313);
plot(t_plot,phi_SS_hat,'r-');grid on;hold on;
ylabel('$$\delta$$$$\phi$$','Interpreter','Latex','FontSize',yfont_size);
plot(t_plot,SM*std_phiSS,'b')
plot(t_plot,-SM*std_phiSS,'b')
ylim([-2*PM3 2*PM3])
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(dEulAx,'xy')
xlim(XLimit)

% True and Estimated Quaternion %
fig=gcf;
figure(fig.Number+1)
QAx(1) = subplot(411);
plot(t_plot,qhat_S_ECI(:,1),'r--');grid on;hold on;
ylabel('$$q_{0}$$','Interpreter','Latex','FontSize',font_size)
plot(t_plot,q_S_ECI(:,1),'b--');
title('Attitude Quaternion','Interpreter','Latex','FontSize',font_size);
legend('Estimate','True','Location','Best')
QAx(2) = subplot(412);
plot(t_plot,qhat_S_ECI(:,2),'r--');grid on;hold on;
ylabel('$$q_{1}$$','Interpreter','Latex','FontSize',font_size)
plot(t_plot,q_S_ECI(:,2),'b--');
QAx(3) = subplot(413);
plot(t_plot,qhat_S_ECI(:,3),'r--');grid on;hold on;
ylabel('$$q_{2}$$','Interpreter','Latex','FontSize',font_size)
plot(t_plot,q_S_ECI(:,3),'b--');
QAx(4) = subplot(414);
plot(t_plot,qhat_S_ECI(:,4),'r--');grid on;hold on;
ylabel('$$q_{3}$$','Interpreter','Latex','FontSize',font_size)
plot(t_plot,q_S_ECI(:,4),'b--');
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(QAx,'xy')
xlim(XLimit)

% Quaternion Estimation Error %
PM1 = 0.05;
PM2 = 0.05;
PM3 = 0.05;
fig=gcf;
figure(fig.Number+1)
dQAx(1) = subplot(311);
plot(t_plot,q_SS_hat(:,2),'r-');grid on;hold on;
ylabel('$$\delta$$$$q_{S\hat{S},1}$$','Interpreter','Latex','FontSize',font_size);
plot(t_plot,SM*std_qSS1,'b-');
plot(t_plot,-SM*std_qSS1,'b-','HandleVisibility','off');
ylim([-2*PM1 2*PM1])
legend('Error',[num2str(SM),'\sigma STD'])
title('Quaternion Estimation Error','Interpreter','Latex','FontSize',font_size);
dQAx(2) = subplot(312);
plot(t_plot,q_SS_hat(:,3),'r-');grid on;hold on;
ylabel('$$\delta$$$$q_{S\hat{S},2}$$','Interpreter','Latex','FontSize',font_size);
plot(t_plot,SM*std_qSS2,'b-');
plot(t_plot,-SM*std_qSS2,'b-','HandleVisibility','off');
ylim([-2*PM2 2*PM2])
dQAx(3) = subplot(313);
plot(t_plot,q_SS_hat(:,4),'r-');grid on;hold on;
ylabel('$$\delta$$$$q_{S\hat{S},3}$$','Interpreter','Latex','FontSize',font_size);
plot(t_plot,SM*std_qSS3,'b-');
plot(t_plot,-SM*std_qSS3,'b-','HandleVisibility','off');
ylim([-2*PM3 2*PM3])
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(dQAx,'xy')
xlim(XLimit)

% Estimated Gyro Biases %
fig=gcf;
figure(fig.Number+1)
BAx(1) = subplot(311);
plot(t_plot,B_true(:,1)*r2d,'b');hold on;
plot(t_plot,b(:,1)*r2d,'r');grid on;
ylabel('$$b_{g,1}$$','Interpreter','Latex','FontSize',font_size);
legend('True Bias','Estimated Bias','Location','Best')
title('Bias Estimates (deg/s)','Interpreter','Latex','FontSize',font_size);
BAx(2) = subplot(312);
plot(t_plot,B_true(:,2)*r2d,'b');hold on;
plot(t_plot,b(:,2)*r2d,'r');grid on;
ylabel('$$b_{g,2}$$','Interpreter','Latex','FontSize',font_size);
BAx(3) = subplot(313);
plot(t_plot,B_true(:,3)*r2d,'b');hold on;
plot(t_plot,b(:,3)*r2d,'r');grid on;
ylabel('$$b_{g,3}$$','Interpreter','Latex','FontSize',font_size);
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(BAx,'xy')
xlim(XLimit)

% Estimated Gyro bias Errors %
PM1 = 0.025;
PM2 = 0.025;
PM3 = 0.025;
fig=gcf;
figure(fig.Number+1)
dQAx(1) = subplot(311);
plot(t_plot,del_b(:,1)*r2d,'r-');grid on;hold on;
ylabel('$$\delta$$$$b_{g,1}$$','Interpreter','Latex','FontSize',font_size);
plot(t_plot,SM*std_db1*r2d,'b-');
plot(t_plot,-SM*std_db1*r2d,'b-','HandleVisibility','off');
ylim([-2*PM1 2*PM1])
legend('Error',[num2str(SM),'\sigma STD'])
title('Bias Estimation Error (deg/s)','Interpreter','Latex','FontSize',font_size);
dQAx(2) = subplot(312);
plot(t_plot,del_b(:,2)*r2d,'r-');grid on;hold on;
ylabel('$$\delta$$$$b_{g,2}$$','Interpreter','Latex','FontSize',font_size);
plot(t_plot,SM*std_db2*r2d,'b-');
plot(t_plot,-SM*std_db2*r2d,'b-','HandleVisibility','off');
ylim([-2*PM2 2*PM2])
dQAx(3) = subplot(313);
plot(t_plot,del_b(:,3)*r2d,'r-');grid on;hold on;
ylabel('$$\delta$$$$b_{g,3}$$','Interpreter','Latex','FontSize',font_size);
plot(t_plot,SM*std_db3*r2d,'b-');
plot(t_plot,-SM*std_db3*r2d,'b-','HandleVisibility','off');
ylim([-2*PM3 2*PM3])
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(dQAx,'xy')
xlim(XLimit)

% Estimated SF %
SFTrue=G_SF;
fig=gcf;
figure(fig.Number+1)
BAx(1) = subplot(311);
plot(t_plot,SFTrue(1)*ones(drl,1),'b');hold on;
plot(t_plot,sf(:,1),'r');grid on;
ylabel('$$sf_{g,1}$$','Interpreter','Latex','FontSize',font_size);
legend('True sf','Estimated sf')
title('Scale Factor Estimates','Interpreter','Latex','FontSize',font_size);
BAx(2) = subplot(312);
plot(t_plot,SFTrue(2)*ones(drl,1),'b');hold on;
plot(t_plot,sf(:,2),'r');grid on;
ylabel('$$sf_{g,2}$$','Interpreter','Latex','FontSize',font_size);
BAx(3) = subplot(313);
plot(t_plot,SFTrue(3)*ones(drl,1),'b');hold on;
plot(t_plot,sf(:,3),'r');grid on;
ylabel('$$sf_{g,3}$$','Interpreter','Latex','FontSize',font_size);
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(BAx,'xy')
xlim(XLimit)

% Estimated SF  Errors %
PM1 = 0.025;
PM2 = 0.025;
PM3 = 0.025;
fig=gcf;
figure(fig.Number+1)
dQAx(1) = subplot(311);
plot(t_plot,del_sf(:,1),'r-');grid on;hold on;
ylabel('$$\delta$$$$sf_{g,1}$$','Interpreter','Latex','FontSize',font_size);
plot(t_plot,SM*std_dsf1,'b-');
plot(t_plot,-SM*std_dsf1,'b-','HandleVisibility','off');
ylim([-2*PM1 2*PM1])
legend('Error',[num2str(SM),'\sigma STD'])
title('Scale Factor Estimation Error','Interpreter','Latex','FontSize',font_size);
dQAx(2) = subplot(312);
plot(t_plot,del_sf(:,2),'r-');grid on;hold on;
ylabel('$$\delta$$$$sf_{g,2}$$','Interpreter','Latex','FontSize',font_size);
plot(t_plot,SM*std_dsf2,'b-');
plot(t_plot,-SM*std_dsf2,'b-','HandleVisibility','off');
ylim([-2*PM2 2*PM2])
dQAx(3) = subplot(313);
plot(t_plot,del_sf(:,3),'r-');grid on;hold on;
ylabel('$$\delta$$$$sf_{g,3}$$','Interpreter','Latex','FontSize',font_size);
plot(t_plot,SM*std_dsf3,'b-');
plot(t_plot,-SM*std_dsf3,'b-','HandleVisibility','off');
ylim([-2*PM3 2*PM3])
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
linkaxes(dQAx,'xy')
xlim(XLimit)

% Angular Error between Estimate and True Spacecraft 3 axis
fig=gcf;
figure(fig.Number+1)
plot(t_plot,Epsilon,'r')
grid on; ylabel('$$\epsilon$$ (deg)','Interpreter','Latex','FontSize',font_size);
hold on
plot(t_plot,std_epsilon*SM,'b')
legend('\epsilon',[num2str(SM),'\sigma STD'])
xlabel(['Time (',timeval,')'],'Interpreter','Latex','FontSize',font_size);
title('Boresight Pointing Error','Interpreter','Latex','FontSize',font_size)
ylim([0 10])
xlim(XLimit)

% Measurement Error Vector(mag)
fig=gcf;
figure(fig.Number+1)
plot(dz_plot(:,1))
hold on; grid on;
plot(dz_plot(:,2))
plot(dz_plot(:,3))
legend('\deltaz_1','\deltaz_2','\deltaz_3')
ylabel('$$\delta$$z (mG)','Interpreter','Latex','FontSize',yfont_size);
xlabel('Measurement','Interpreter','Latex','FontSize',yfont_size)
title('Measurement Error Vector, TAM','Interpreter','Latex','FontSize',font_size)
xlim(XLimitdz)

% Measurement Error Vector (sun)
if const.UseSun ==1
    fig=gcf;
    figure(fig.Number+1)
    plot(dz_plot(:,4))
    hold on; grid on;
    plot(dz_plot(:,5))
    plot(dz_plot(:,6))
    legend('\deltaz_1','\deltaz_2','\deltaz_3')
    ylabel('$$\delta$$z','Interpreter','Latex','FontSize',yfont_size);
    xlabel('Measurement','Interpreter','Latex','FontSize',yfont_size)
    title('Measurement Error Vector, Sun Sensor','Interpreter','Latex','FontSize',font_size)
    xlim(XLimitdz)
end

if const.RunSOTest ==1
    % Deterministic Observability % 
    fig=gcf;
    figure(fig.Number+1)
    scatter(t_mag,Obs_G_rank,'r*')
    grid on; hold on;
    scatter(t_mag,Obs_G_rank_true,'k.')
    ylabel('Rank \it $$O_{k}^{T}O_{k}$$','Interpreter','Latex','FontSize',font_size);
    xlabel(['Time (sec)'],'Interpreter','Latex','FontSize',font_size);
    legend('Noisy System','True System')
    title('Deterministic Observability','Interpreter','Latex','FontSize',font_size)
    ylim([0 FullRank+2])
    xlim([0 10/const.MagSamp])

    % Stochastic Observability %
    fig=gcf;
    figure(fig.Number+1)
    subplot(211);
    plot(t_mag,L_sv,'r-');grid on; hold on;
    plot(t_mag,L_sv_true,'k--');
    ylabel('$$\sigma_{max}$$($$\Lambda_{k})$$','Interpreter','Latex','FontSize',font_size);
    xlim([0 10*1/const.MagSamp])
    xlabel(['Time (sec)'],'Interpreter','Latex','FontSize',font_size);
    legend('Noisy System','True System')
    title('Stochastic Observability','Interpreter','Latex','FontSize',font_size)
    subplot(212);
    plot(t_mag/3600,Q_sv,'r-');grid on; hold on;
    plot(t_mag/3600,Q_sv_true,'k--');
    ylabel('$$\sigma_{max}$$($$\overline {Q}_{k})$$','Interpreter','Latex','FontSize',font_size);
    ylim([0 5e-3])
    xlabel(['Time (hrs)'],'Interpreter','Latex','FontSize',font_size);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%