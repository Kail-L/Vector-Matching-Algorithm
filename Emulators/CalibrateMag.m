%% Calibrate Mag %%
% This is a script to Batch calibrate magnetometer measurements emulated from the
% emulateMag.m function.
%
% Kail Laughlin
% 3/12/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   A model that relates the measured EMF (EMF_m) to the true EMF
%   (EMF_true) is (note: the scalar components of the vector magnetometer
%   calibration equation are expressed in the same frame)

%   EMF_m = C_sisfm*EMF_true + b_hardironbias +b_whitenoise

%       where

%   C_sisfm = C_Misalignmenterror*C_ScaleFactorerror*C_SoftIronerror

%   C_sisfm is a 3x3 matrix and is the product of three 3x3 matrices

EMF_m_gbl = [];
EMF_cal_gbl = [];

EMF_m_check = [];
EMF_cal_check = [];

emf_count = 1;
len_Mag_time = length(t_mag);
for i = 2:len_Mag_time
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % select the measured EMF data for the current SGPB position and time:
    emf_meas = BbNoisy(i,:)';
    % select the EMF model data for the current SGPB position and time:
    emf_model = Bb(i,:)';
    % formulate the time-to-date EMF measurement matrix:
    EMF_m_gbl = [EMF_m_gbl;
                 emf_meas];
    % formulate the EMF measurement matrix for checking the time dependence
    % of the calibration parameters
    EMF_m_check = [EMF_m_check;
                   emf_meas];
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % formulate the EMF calibration matrix for the current SGPB position
    % and time:
    
    emf_cal_mat = [emf_model' zeros(1,3) zeros(1,3) 1 0 0;
                   zeros(1,3) emf_model' zeros(1,3) 0 1 0;
                   zeros(1,3) zeros(1,3) emf_model' 0 0 1];
    
    % formulate the time-to-date EMF calibration matrix
    
    EMF_cal_gbl = [EMF_cal_gbl;
                   emf_cal_mat];
    
    % formulate the EMF calibration matrix for checking the time dependence
    % of the calibration parameters
    
    EMF_cal_check = [EMF_cal_check;
                     emf_cal_mat]; 
    if (emf_count == 5000)
        % compute the calibration parameters using the current time-to-date
        % calibration matrix
        
        emf_count
        EMF_cal_par_check = inv(EMF_cal_check'*EMF_cal_check)*...
            EMF_cal_check'*EMF_m_check;
                
        % reset the check parameters
        
        EMF_m_check = [];
        EMF_cal_check = [];
        emf_count = 0;
    end
    emf_count = emf_count + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determine the calibration parameters using the global EMF calibration 
% matrix
    
emf_cal_par = inv(EMF_cal_gbl'*EMF_cal_gbl)*EMF_cal_gbl'*EMF_m_gbl;
EMF_cal_par_gbl = emf_cal_par;

% determine the calibrated magnetic field measurements using the current
% EMF calibration parameters:
    
emf_sisfm = [emf_cal_par(1) emf_cal_par(2) emf_cal_par(3);
             emf_cal_par(4) emf_cal_par(5) emf_cal_par(6);
             emf_cal_par(7) emf_cal_par(8) emf_cal_par(9)];
emf_sisfm_inv = inv(emf_sisfm);

emf_bias = [emf_cal_par(10);
            emf_cal_par(11);
            emf_cal_par(12)];

for i = 1:len_Mag_time
    emf_m_cal_gbl = emf_sisfm_inv*(BbNoisy(i,:)' - emf_bias);
    EMF_m_cal(:,i) = emf_m_cal_gbl;
end
EMF_m_cal = EMF_m_cal';

fig=gcf;
figure(fig.Number+1)
plot3(Bb(:,1),Bb(:,2),Bb(:,3),'bo')
grid;ylabel('B (T)');title('Error Free B, Body Frame')

fig=gcf;
figure(fig.Number+1)
plot3(BbNoisy(:,1),BbNoisy(:,2),BbNoisy(:,3),'go')
grid;ylabel('B (T)');title('Corrupted B, Body Frame')

fig=gcf;
figure(fig.Number+1)
plot3(EMF_m_cal(:,1),EMF_m_cal(:,2),EMF_m_cal(:,3),'ro')
grid;ylabel('B (T)');title('Calibrated B, Body Frame')