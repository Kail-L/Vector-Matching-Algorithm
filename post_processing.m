%% Post processing %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to post process data for a satellite ran through the
% SatAttEst.m code. 

% Kail Laughlin
% 7/6/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% r2d=180/pi;
% d2r=pi/180;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract attitude states %%
x_att_out = x_out(:,7:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Get Sun Info %%
% If using sim data:
%  - Get Sun Vector (ECI) using planetEphemeris.m (Matlab function) 
% If using SGPB data:
%  - Get Sun Vector (ECI) from r_S_ECI_uv data in SGPB_Sun_Pos_Data.m

if const.UseExternal==1
    if const.UseSun==1
        S_ECI=ExtVec(:,7:9);
    else
        S_ECI=zeros(length(t_sun),3);
    end
else
    if const.UseSun==1
        S_ECI=planetEphemeris(juliandate(t_utc_sun),'Earth','Sun');
    else 
        S_ECI=zeros(length(t_sun),3);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Matrices for s p e e d %%
B_ECI=zeros(length(t_mag),3);      % Mag Field Vector (ECI)
B_S=zeros(length(t_mag),3);      % Mag Field Vector, clean (Body)
B_S_noisy=zeros(length(t_mag),3);     % Calibrated mag field vector (only for SGPB data)
S_SpaceECI=zeros(length(t_sun),3); % Sun Vector to Spacecraft, clean (Body)
S_S=zeros(length(t_sun),3); % Sun Vector, clean (Body)
S_S_noisy=zeros(length(t_sun),3);  % Sun Vector, noisy (Body)
R_ECI=zeros(length(t),3);      % Position Vector (ECI)
R_S=zeros(length(t),3);      % Position Vector (Body)
V_ECI=zeros(length(t),3);      % Velocity Vector (ECI)
V_S=zeros(length(t),3);      % Position Vector (Body)
q_S_ECI_true = zeros(length(t),4);  % True quaternion from inertial to body frame
norm_q = zeros(length(t),1);% Norm of true quaternion - 1 for check
Euler_angles=zeros(length(t),3);    % Euler Angles of SC
Euler_angles_rg=zeros(length(t),3);    % Euler Angles of SC
ddotEuler_angles=zeros(length(t),3);% Rate of change of Euler Angles
omega_S_ECI_true=zeros(length(t),3);        % p,q,r of SC
RateM=zeros(3,3);                   % Rate matrix used to calc deriv of Euler Angles
E = zeros(length(t),1);             % Total Energy of system
magcount=1;
suncount=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run Main Post Processing Loop %%
for lv1 = 1:length(x_out)
    r_ECI = x_out(lv1,1:3)';   % Position
    R_ECI(lv1,:)=r_ECI;           % Store Position Vector (ECI)
    v_ECI = x_out(lv1,4:6)';   % Velocity
    V_ECI(lv1,:)=v_ECI;           % Store Velocity Vector (ECI)
    r_norm = sqrt(r_ECI'*r_ECI); % Norm of position

    q_S_ECI_true(lv1,:) = x_att_out(lv1,1:4); % True quaternion

    omega_S_ECI_true(lv1,:) = x_att_out(lv1,5:7); % Angular velocity of b frame axes, resolved in b frame
    CS_ECI = quat2dcm(q_S_ECI_true(lv1,:));          % DCM from inertial a frame to body b frame
    eul_S_ECI = dcm2eul(CS_ECI);                % Euler angles 
    psi_S_ECI = eul_S_ECI(1);
    the_S_ECI = eul_S_ECI(2);
    phi_S_ECI = eul_S_ECI(3);

    % Determine rate matrix for this instant in time %
    RateM=[1 sin(phi_S_ECI)*tan(phi_S_ECI) cos(phi_S_ECI)*tan(phi_S_ECI);
           0 cos(phi_S_ECI) -sin(phi_S_ECI);
           0 sin(phi_S_ECI)/cos(the_S_ECI) cos(phi_S_ECI)/cos(the_S_ECI)];

    % Store Euler Angle and Euler Rate values %
    ddotEulAng = RateM*omega_S_ECI_true(lv1,:)';          % Euler Angle rates of change
    Euler_angles(lv1,:) = [phi_S_ECI the_S_ECI psi_S_ECI]; % Store Euler angles
    ddotEuler_angles(lv1,:)= [ddotEulAng(1),ddotEulAng(2),ddotEulAng(3)]; 

    % Check Sim Through Quaternion Norm calculations %
    norm_q(lv1) = q_S_ECI_true(lv1,:)*q_S_ECI_true(lv1,:)' - 1; % Norm of quaternion minus 1

    % Frame Transformations %
    r_S = CS_ECI*r_ECI;   % Position in body frame
    v_S = CS_ECI*v_ECI;   % Vel in body frame
    R_S(lv1,:) = r_S; % Store position vector
    V_S(lv1,:) = v_S; % Store Vel Vector  

    % Mag Field %
    if const.UseMeas && magcount<=length(t_mag)&&t(lv1)==t_mag(magcount)
        tBA = t_utc(lv1,:);
        ST = datevec(SimStart);
        B_ECI(magcount,:) = GetMagData(ST(1:3),r_ECI,tBA,t_mag(magcount),const);
        B_S(magcount,:)= CS_ECI*B_ECI(magcount,:)';
        magcount=magcount+1;
    end

    % Sun Vector %
    if const.UseSun && suncount<=length(t_sun)&&t(lv1)==t_sun(suncount)
        if const.UseSun==1
            S_SpaceECI(suncount,:)=S_ECI(suncount,:)-(r_ECI'./1000);
            S_ECI(suncount,:) = S_ECI(suncount,:)/norm(S_ECI(suncount,:)); % Normalize sun Meas.
            theta_sun = acos(dot(S_ECI(suncount,:),r_ECI)/(norm(S_ECI(suncount,:))*r_norm));
            phi_sun = asin(const.Re/r_norm); % Eclipsing angle 
            S_ECI(suncount,:) = S_SpaceECI(suncount,:)/norm(S_SpaceECI(suncount,:));
            if theta_sun<=phi_sun
                S_ECI(suncount,:) = [0 0 0];
            end 
            S_S(suncount,:)=CS_ECI*S_ECI(suncount,:)';
            [S_S_noisy(suncount,:),~]=SunSensorNoisy(S_ECI(suncount,:),CS_ECI,const);
        else
            S_S_noisy(suncount,:) = [0 0 0];
        end
        suncount=suncount+1;
    end

    % Total Energy with J2 perturbation %
    E(lv1) = 0.5*const.m_s*(v_ECI'*v_ECI) - const.mu1*const.m_s/r_norm...
        + const.mu1*const.m_s/r_norm^3*const.J2*const.Re^2*(3/2*(r_ECI(3)/r_norm)^2-0.5);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Gyro Parameters %%
Gyro_Good = [omega_S_ECI_true(:,1) omega_S_ECI_true(:,2) omega_S_ECI_true(:,3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Mag Parameters %%
if const.UseSGPB ==1
    B_ECI = ExtVec(:,1:3);
    B_S = ExtVec(:,4:6);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rename Parameters %%
phi=Euler_angles(:,1);
theta=Euler_angles(:,2);
psi=Euler_angles(:,3);
phi_rg=Euler_angles_rg(:,1);
theta_rg=Euler_angles_rg(:,2);
psi_rg=Euler_angles_rg(:,3);
phidot=ddotEuler_angles(:,1);
thetadot=ddotEuler_angles(:,2);
psidot=ddotEuler_angles(:,3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Package Data %%
true_att=[t psi psidot theta thetadot phi phidot Gyro_Good];

%% Create Vectors of Data %%
Pos=[t R_ECI R_S];
Vel=[t V_ECI V_S];
Mag=[B_ECI B_S]; 
Sun=[S_ECI S_S S_S_noisy]; % Unit Vectors
IMU=[t, Gyro_Good];
TrueAng=[t psi psidot theta thetadot phi phidot];
RangeAng=[t psi_rg theta_rg phi_rg];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% Save Data %%
% if const.UseSun==1
%     save('IMU_Eul_SV_Mag.mat','IMU','TrueAng','RangeAng','Mag','Sun')
% else
%     save('IMU_Eul_SV_Mag.mat','IMU','TrueAng','RangeAng','Mag')
% end
% 
% save('SatPosAtt.mat','true_att','Pos','Vel')
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
