%% Stochastic Observability Test for VMA_EKF %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to determine stochastic observability (as defined in ch.3 of
% Vibhor Bageshwar's 2008 Ph.D Thesis) of small satellite (cubesat). Script
% takes EKF data from VMA_EKF.m and performs the stochastic observability
% algorithm on it. 
%
% Note: Script calculates stochastic observability for two systems: the
% system created by linearizing about the truth trajectory, and the system
% created by linearizing about the estimate trajectory.
%
% See ch. 3 of Vibhor Bageshwar's 2008 Thesis for detailed information on
% stochastic observability algorithm.
%
% Kail Laughlin, Vibhor Bageshwar
% 7/8/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define Constants %%
piTol = 1e-3;    % pseudo inverse tolerance
piFlag = 0;       % 0: don't use tolerance, 1: use tolerance
PHI_gbl = eye(PHI_row,PHI_col); % Global state transition matrix
PHI_c = eye(PHI_row,PHI_col);   % Intermediate state transition matrix
PHI_int = eye(PHI_row,PHI_col); % State transition matrix between measurements
Qd_c = zeros(PHI_row,PHI_col);
W_l = eye(PHI_row,PHI_col);     % null space projections
W_r = eye(PHI_row,PHI_col);     % null space projections
OBS_mat = [];                   % deterministic obsv. matrix
OBS_gram_test = zeros(PHI_row,PHI_col);
Tol_rank = 1e-10;               % SVD tolerance for rank computations
mu_counter = 2;                 % measurement update counter
Iter_num = TU_num;              % number of iterations to run
IC_flag = 0;                    % flag to determine if P depends on initial conditions
SOFlag = 0;                     % flag used to state if Q_bar grows infinitely
tu_per_mu = const.SampRatio;    % time update to measurement update ratio
fig=gcf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize Matrices for s p e e d %%
Lambda_sv = zeros(PHI_col,MU_num);    % singular values of Lambda
Q_bar_sv = zeros(PHI_col,MU_num);     % singular values of Q_bar
Lambda_std = zeros(PHI_col,MU_num);   % standard deviations of diagonals of Lambda 
Q_bar_std = zeros(PHI_col,MU_num);    % standard deviations of Q_bar
r_OBS_gram = zeros(1,MU_num);         % rank of observability grammian
OBS_gram_sv = zeros(PHI_col,MU_num);  % singular values of observability grammian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Run SO Test %%
for i=2:Iter_num
%%% Time update: Get matrices for current time %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if tu_per_mu == 1
        PHI_c = PHI_SO{i};          % Current discreet F matrix
        Qd_c = Qd_SO{i};            % Current discreet PSD matrix
        PHI_int = PHI_c;
    elseif tu_per_mu > 1 % time update rate > measurement update rate
        Qd_c = Qd_SO{i} + PHI_int*Qd_c*PHI_int';
        PHI_c = PHI_SO{i};
        PHI_int = PHI_c*PHI_int;
    end
    
%%% Measurement update %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if const.UseMeas && mu_counter <= MU_num && TU_time(i) == MU_time(mu_counter)  
        
        R = Rz_SO{mu_counter};  % Measurement noise covariance matrix
        H = H_SO{mu_counter};   % Current Measurement matrix
        
%%%%%%% SO Test %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if IC_flag == 0         % P depends on initial conditions
            if mu_counter == 2  % First measurement update    
                [RU,RS,RV]=svd(R);
                Gamma = RU*(sqrt(RS))*RU';  % Matrix Square Root of Gamma
                [GU,GS,GV]=svd(Gamma); 
                rankGinvout = 0;
                iGS = zeros(H_row,H_row);
                iGamma = piTol*eye(H_row,H_row);
                for iInv = 1:H_row
                    if GS(iInv,iInv)>piTol
                        rankGinvout = iInv;
                        iGS(iInv,iInv) = 1/GS(iInv,iInv);
                    end                  
                end
                iGamma = GV(:,1:rankGinvout)*iGS(1:rankGinvout,1:rankGinvout)*GV(:,1:rankGinvout)';
           
                Theta = H'*iGamma;
                ThT_Th = Theta'*Theta;
                
                % Calculate pseudo inverse of Theta^T*Theta %
%                 [ThU,ThS,ThV]=svd(ThT_Th);      % Matrix Pseudoinverse
%                 rankThinvout = 0;
%                 iThS = zeros(H_row,H_row);
%                 Theta_pi = ThS;%piTol*eye(H_row,H_row);
%                 for iInv = 1:H_row
%                     if ThS(iInv,iInv)>piTol
%                         rankThinvout = iInv;
%                         iTh(iInv,iInv) = i/ThS(iInv,iInv);
%                     end
%                 end
%                 Theta_pi = ThV(:,1:rankThinvout)*iTh(1:rankThinvout,1:rankThinvout)*ThU(:,1:rankThinvout)';                
                [ThU,ThS,ThV]=svd(ThT_Th);
                for i = 1:H_row
                    for j=1:H_row
                        if abs(ThS(i,j))<piTol
                            iThS(i,j) = 0;                            
                        else
                            iThS(i,j) = 1/ThS(i,j);
                        end  
                    end
                end                

                Theta_pi = ThV*iThS*ThV';
                
                % Omega_centered value from nullspace projection %
                W_c = eye(PHI_row,PHI_col) - Theta*Theta_pi*Theta';
                
                % Lambda %
                Lambda = PHI_c*PHI_gbl*W_c*PHI_gbl'*PHI_c';
                Lambda_sv(:,mu_counter) = svd(Lambda); % singular values
                
                for k = 1:PHI_row % Make singular values zero if below a threshold
                    X = Lambda;
                    if abs(X(k,k))<0
                        X(k,k) = 0;
                    end
                    Lambda_std(k,mu_counter) = sqrt(abs(X(k,k))); % std deviation
                end
                
                % Q bar %
                Q_bar = Qd_c + PHI_c*PHI_gbl*Theta*(Theta_pi*Theta_pi)*Theta'*PHI_gbl'*PHI_c';  
                Q_bar_sv(:,mu_counter) = svd(Q_bar); % singular values
                
                for k = 1:PHI_row 
                    X = Q_bar;
                    if abs(X(k,k))<0
                        X(k,k) = 0;
                    end
                    Q_bar_std(k,mu_counter) = sqrt(abs(X(k,k))); % std deviation
                end
                
                % Deterministic observability
                OBS_mat = H;
                OBS_gram = OBS_mat'*OBS_mat;
                c_OBS_gram(mu_counter) = cond(OBS_gram);
                r_OBS_gram(mu_counter) = rank(OBS_gram);
                OBS_gram_sv(:,mu_counter) = svd(OBS_gram);
                
            else % remaining measurement updates
                TG = sqrtm(R + H*Q_bar*H');
                [RU,RS,RV]=svd(R + H*Q_bar*H');
                Gamma = RU*(sqrt(RS))*RU';  % Matrix Square Root of Gamma
                
                [GU,GS,GV]=svd(Gamma);      % Matrix Pseudoinverse
                rankGinvout = 0;
                iGS = zeros(H_row,H_row);
                iGamma = piTol*eye(H_row,H_row);
                for iInv = 1:H_row
                    if GS(iInv,iInv)>piTol
                        rankGinvout = iInv;
                        iGS(iInv,iInv) = 1/GS(iInv,iInv);
                    end
                end
                iGamma = GV(:,1:rankGinvout)*iGS(1:rankGinvout,1:rankGinvout)*GV(:,1:rankGinvout)';                      
                Theta = W_r'*PHI_gbl'*H'*iGamma;
                ThT_Th = Theta'*Theta;
                
                % Calculate pseudo inverse of Theta^T*Theta %
                [ThU,ThS,ThV]=svd(ThT_Th);      % Matrix Pseudoinverse
%                 rankThinvout = 0;
%                 iThS = zeros(H_row,H_row);
%                 Theta_pi = piTol*eye(H_row,H_row);
%                 for iInv = 1:H_row
%                     if ThS(iInv,iInv)>piTol
%                         rankThinvout = iInv;
%                         iTh(iInv,iInv) = i/ThS(iInv,iInv);
%                     end
%                 end
%                 Theta_pi = ThV(:,1:rankThinvout)*iTh(1:rankThinvout,1:rankThinvout)*ThV(:,1:rankThinvout)';  
                for i = 1:H_row
                    for j=1:H_row
                        if abs(ThS(i,j))<piTol
                            iThS(i,j) = 0;                            
                        else
                            iThS(i,j) = 1/ThS(i,j);
                        end                 
                    end
                end 
                    
                Theta_pi = ThV*iThS*ThV';
                
                % Omega_centered value from nullspace projection %
                W_c = eye(PHI_row,PHI_col) - Theta*Theta_pi*Theta';
                
                % Lambda %
                Lambda = PHI_c*PHI_gbl*W_l*W_c*W_r*PHI_gbl'*PHI_c';
                Lambda_sv(:,mu_counter) = svd(Lambda); % singular values
                
                for k = 1:PHI_row % Make singular values zero if negative
                    X = Lambda;
                    if abs(X(k,k))<0
                        X(k,k) = 0;
                    end
                    Lambda_std(k,mu_counter) = sqrt(abs(X(k,k))); % std deviation
                end
                
                % Q_bar %

                for i = 1:H_row
                    for j=1:H_row
                        if abs(ThT_Th(i,j))<0.001
                            ThT_Th(i,j) = 0;
                        end
                    end
                end
                
                % Calculate pseudo inverse of Theta^T*Theta %
                [ThU,ThS,ThV]=svd(ThT_Th);      % Matrix Pseudoinverse
%                 rankThinvout = 0;
%                 iThS = zeros(H_row,H_row);
%                 Theta_pi = piTol*eye(H_row,H_row);
%                 for iInv = 1:H_row
%                     if ThS(iInv,iInv)>piTol
%                         rankThinvout = iInv;
%                         iTh(iInv,iInv) = i/ThS(iInv,iInv);
%                     end
%                 end
%                 Theta_pi = ThV(:,1:rankThinvout)*iTh(1:rankThinvout,1:rankThinvout)*ThV(:,1:rankThinvout)';                 
                for i = 1:H_row
                    for j=1:H_row
                        if abs(ThS(i,j))<piTol
                            iThS(i,j) = 0;                            
                        else
                            iThS(i,j) = 1/ThS(i,j);
                        end                 
                    end
                end     
                Theta_pi = ThV*iThS*ThV';
                 
                % Determine nullspace projector? Ask Vibhor %
                [U,Sig,V] = svd(Theta);
                r_SigT_Sig = rank(Sig'*Sig);  %rank (sigma^T sigma)
                [SigT_Sig_row,SigT_Sig_col] = size(Sig'*Sig); % number of rows and columns
                proj_SigT_Sig = zeros(SigT_Sig_row,SigT_Sig_col); % set projector to zero
                if SigT_Sig_row-r_SigT_Sig>0
                    proj_SigT_Sig(r_SigT_Sig+1:SigT_Sig_row,r_SigT_Sig+1:SigT_Sig_col) = ...
                        eye(SigT_Sig_row-r_SigT_Sig,SigT_Sig_col-r_SigT_Sig);
                end
                
                Q_bar = PHI_int*Q_bar*PHI_int' + Qd_c...
                    + PHI_c*PHI_gbl*W_l*Theta*(Theta_pi*Theta_pi)*Theta'*W_r*PHI_gbl'*PHI_c' ...
                    - PHI_c*PHI_gbl*W_l*Theta*Theta_pi*iGamma'*H*Q_bar*PHI_int' ...
                    - PHI_int*Q_bar*H'*iGamma*Theta_pi*Theta'*W_r*PHI_gbl'*PHI_c' ...
                    - PHI_int*Q_bar*H'*iGamma*V*proj_SigT_Sig*V'*iGamma'*H*Q_bar*PHI_int';
                               
                QStor{mu_counter}=Q_bar;
                Q_bar_sv(:,mu_counter) = svd(Q_bar); % singular values  
                
                for k = 1:PHI_row 
                    X = Q_bar;
                    if abs(X(k,k))<0
                        X(k,k) = 0;
                    end
                    Q_bar_std(k,mu_counter) = sqrt(abs(X(k,k))); % std deviation
                end
                
                % deterministic observability %
                OBS_mat = [OBS_mat;
                           H*PHI_gbl];
                OBS_gram = OBS_mat'*OBS_mat;
                c_OBS_gram(mu_counter) = cond(OBS_gram);
                r_OBS_gram(mu_counter) = rank(OBS_gram);
                OBS_gram_sv(:,mu_counter) = svd(OBS_gram);
            end
            
            % check if P still a function of initial conditions
            if Lambda_sv(1,mu_counter) < Tol_rank % P is not a function of initial conditions
                IC_flag = 1;
            end
                
            % Global matrices
            W_l = W_l*W_c; % Omega_0 Omega_1 ... Omega_k
            W_r = W_c*W_r; % Omega_k ... Omega_1 Omega_0
            
        else % P does not depend on initial conditions

            [RU,RS,RV]=svd(R + H*Q_bar*H');
            Gamma = RU*(sqrt(RS))*RU';  % Matrix Square Root of Gamma
            [GU,GS,GV]=svd(Gamma); 
            rankGinvout = 0;
            iGS = zeros(H_row,H_row);
            iGamma = piTol*eye(H_row,H_row);
            for iInv = 1:H_row
                if GS(iInv,iInv)>piTol
                    rankGinvout = iInv;
                    iGS(iInv,iInv) = 1/GS(iInv,iInv);
                end
            end
            iGamma = GV(:,1:rankGinvout)*iGS(1:rankGinvout,1:rankGinvout)*GV(:,1:rankGinvout)';            
%             for i = 1:H_row
%                 for j=1:H_row
%                     if abs(GS(i,j))<piTol
%                         iGS(i,j) = 0;                            
%                     else
%                         iGS(i,j) = 1/GS(i,j);
%                     end                 
%                 end
%             end                 
%             iGamma = GV*iGS*GU';             

            % Q_bar %

            Q_bar = PHI_int*Q_bar*PHI_int' + Qd_c...
                - PHI_int*Q_bar*H'*iGamma*iGamma'*H*Q_bar*PHI_int';

            Q_bar_sv(:,mu_counter) = svd(Q_bar);  
            QStor{mu_counter}=Q_bar;

            for k = 1:PHI_row 
                X = Q_bar;
                if abs(X(k,k))<0
                    X(k,k) = 0;
                end
                Q_bar_std(k,mu_counter) = sqrt(abs(X(k,k)));

            end

            if Q_bar_sv(1,mu_counter)>1e8
                SOFlag=1;
            end

            Q_bar_sv(:,mu_counter) = svd(Q_bar); 


            % deterministic observability
            OBS_mat = [OBS_mat;
                       H*PHI_gbl];
            OBS_gram = OBS_mat'*OBS_mat;
            c_OBS_gram(mu_counter) = cond(OBS_gram);
            r_OBS_gram(mu_counter) = rank(OBS_gram);
            OBS_gram_sv(:,mu_counter) = svd(OBS_gram);   
        end
        
        if tu_per_mu == 1
            PHI_gbl = PHI_c*PHI_gbl;
        end
        mu_counter = mu_counter+1;
        Qd_c = zeros(PHI_row,PHI_col); % Reset Qd_c value after measurement
        PHI_int = eye(PHI_row,PHI_col);% Reset state transition matrix between measurements
    end
    if tu_per_mu > 1
        PHI_gbl = PHI_c*PHI_gbl;      
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                             