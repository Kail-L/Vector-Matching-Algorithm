function [M_cal, D_new, b_bold_new, P_new, thetaP_new] = CalibrateTAM(M_noisy, M_model, thetaPrime, Ptheta)
%% Magnetometer Calibration %%
% Function to calibrate Magnetometer data output from emulateMag.m.

% Measurement Covariance %
% Scalar value, representing model accuracy.

R_mag = 0.5; % milliGauss

% State Vector %
c = thetaPrime(1:3);
E_bold = thetaPrime(4:end);
E = [E_bold(1) E_bold(4) E_bold(5);
     E_bold(4) E_bold(2) E_bold(6);
     E_bold(5) E_bold(6) E_bold(3)];
 
[U, V] = eig(E);
W = diag([-1+sqrt(1+V(1,1)) -1+sqrt(1+V(2,2)) -1+sqrt(1+V(3,3))]);
D = U*W*U';
b_bold = inv(eye(3)+D)*c;

B_k = M_noisy';
B_mod_k = M_model'; %model in ECI frame
S_k = [B_k(1)^2 B_k(2)^2 B_k(3)^2 2*B_k(1)*B_k(2) 2*B_k(1)*B_k(3) 2*B_k(2)*B_k(3)];
L_k = [2*B_k' -S_k];
h_k = L_k*thetaPrime - norm(b_bold)^2;
Ec = inv(eye(3) + E)*c;
dbdc = 2*Ec;
for m =1:3
    for n=1:3
        if m == n
            Kdelta = 1;
        else
            Kdelta = 0;
        end
        dbde(m,n) = -(2-Kdelta)*Ec(m)*Ec(n);
    end
end
dbdE = [dbde(1,1) dbde(2,2) dbde(3,3) dbde(1,2) dbde(1,3) dbde(2,3)];
dbdThetaP = [dbdc' dbdE];
H_cal = L_k - dbdThetaP;
zk = norm(B_k)^2-norm(B_mod_k)^2;
K = (Ptheta*H_cal')/(H_cal*Ptheta*H_cal' + R_mag);
x = thetaPrime + K*(zk - h_k);
P = (eye(9) - K*H_cal)*Ptheta;

% Recompute bias and D matrix %
thetaP_new = x;
c_new = x(1:3);
E_new = x(4:end);
E = [E_new(1) E_new(4) E_new(5);
     E_new(4) E_new(2) E_new(6);
     E_new(5) E_new(6) E_new(3)];
[U, V] = eig(E);
W = diag([-1+sqrt(1+V(1)) -1+sqrt(1+V(2)) -1+sqrt(1+V(3))]);
D_new = U*W*U';
b_bold_new = inv(eye(3)+D_new)*c_new;
McD = [b_bold_new(1) 0 0 b_bold_new(2) b_bold_new(3) 0;
       0 b_bold_new(2) 0 b_bold_new(1) 0 b_bold_new(3);
       0 0 b_bold_new(3) 0 b_bold_new(1) b_bold_new(2)];
MED = [2*D_new(1,1) 0 0 2*D_new(1,2) 2*D_new(1,3) 0;
       0 2*D_new(2,2) 0 2*D_new(1,2) 0 2*D_new(2,3);
       0 0 2*D_new(3,3) 0 2*D_new(1,3) 2*D_new(2,3);
       D_new(1,2) D_new(1,2) 0 D_new(1,1)+D_new(2,2) D_new(2,3) D_new(1,3);
       D_new(1,3) 0 D_new(1,3) D_new(2,3) D_new(1,1)+D_new(3,3) D_new(1,3);
       0 D_new(2,3) D_new(2,3) D_new(1,3) D_new(1,2) D_new(2,2)+D_new(3,3)];
JacobianDE = inv([(eye(3) + D_new) McD;
              zeros(6,3) 2*eye(6) + MED]);
P_new = JacobianDE*P*JacobianDE';

M_cal = (eye(3)+D_new)*M_noisy' - b_bold_new;

end
