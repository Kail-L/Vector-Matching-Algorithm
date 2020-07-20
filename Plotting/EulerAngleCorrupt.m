%% EulerAngleCorrupt
% Script to get corrupt and true euler angle plots 
EulAng0 = [0 0 0]';         % Initial True Euler Angles
EulAngGood(1,:)=EulAng0;
EulAngBad(1,:)=EulAng0;

for i=2:length(t)
    PsiGood=EulAngGood(i-1,1);
    ThetaGood=EulAngGood(i-1,2);
    PhiGood=EulAngGood(i-1,3);
    RateGood=[1 sin(PhiGood)*tan(PhiGood) cos(PhiGood)*tan(PhiGood);
       0 cos(PhiGood) -sin(PhiGood);
       0 sin(PhiGood)/cos(ThetaGood) cos(PhiGood)/cos(ThetaGood)];
    dotEAngGood=RateGood*Gyro_Good(i-1,:)';
   
    PsiBad=EulAngBad(i-1,1);
    ThetaBad=EulAngBad(i-1,2);
    PhiBad=EulAngBad(i-1,3);
    RateBad=[1 sin(PhiBad)*tan(PhiBad) cos(PhiBad)*tan(PhiBad);
       0 cos(PhiBad) -sin(PhiBad);
       0 sin(PhiBad)/cos(ThetaBad) cos(PhiBad)/cos(ThetaBad)];
    dotEAngBad=RateBad*Gyro_Corrupt(i-1,:)';
    
    EulAngGood(i,:)=EulAngGood(i-1,:)+(dotEAngGood*dt)';
    EulAngBad(i,:)=EulAngBad((i-1),:)+(dotEAngBad*dt)';
end

% Plot Omega Values
fig=gcf;
figure(fig.Number+1)
subplot(311)
plot(t,Gyro_Corrupt(:,1),'r-')
hold on
grid on
plot(t,Gyro_Good(:,1),'k-')
ylabel('\omega_S_1 (rad/s)','FontSize',font_size)
legend('Noisy Gyro','Clean Gyro')
subplot(312)
plot(t,Gyro_Corrupt(:,2),'r-')
hold on
grid on
plot(t,Gyro_Good(:,2),'k-')
ylabel('\omega_S_2 (rad/s)','FontSize',font_size)
subplot(313)
plot(t,Gyro_Corrupt(:,3),'r-')
hold on
grid on
plot(t,Gyro_Good(:,3),'k-')
ylabel('\omega_S_3 (rad/s)','FontSize',font_size)
xlabel('Time (sec)','FontSize',font_size)

% Plot psi
fig=gcf;
figure(fig.Number+1)
plot(thr,yaw_hat_ba,'r-',thr,yaw_ba,'k-');grid on;ylabel('\psi (deg)','FontSize',font_size);
xlabel('Time (sec)','FontSize',font_size)
title('Yaw Propagation w/ Clean & Noisy Gyro')
legend('Noisy Gyro','Clean Gyro')

   