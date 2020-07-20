%% Test Script to run calibration function %%

uECI=Mag(:,1:3)*T2mg;          % Mag field vector in inertial frame, mGauss
Mag_Clean_S=Mag(:,4:6)*T2mg;   % Clean mag in body frame

[uS, CTrue, Btrue] = emulateMag(t_mag,Mag_Clean_S,const); % Noisy mag measurement

tsmag = length(t_mag);
P0 = [500*eye(3) zeros(3,6);
      zeros(6,3) 0.001*eye(6)];
P = P0;
thetaP = zeros(9,1);
for i=1:tsmag
    [M_cal, D_new, b_bold_new, P_new, thetaP_new] = CalibrateTAM(uS(i,:), uECI(i,:), thetaP, P);
    CalMag(i,:) = M_cal;
    Best(i,:) = b_bold_new;
    Cest{i} = inv(eye(3)+D_new);
    Pstor{i} = P_new;
    thetaP = thetaP_new;
    P = P_new;
end

% Plot Results %
subplot(311)
plot(t_mag,Btrue(:,1),'b');hold on; grid on;
plot(t_mag,Best(:,1),'r');
legend('True Hard Iron Bias','Estimated Hard Iron Bias')
subplot(312)
plot(t_mag,Btrue(:,2),'b');hold on; grid on;
plot(t_mag,Best(:,2),'r');
subplot(313)
plot(t_mag,Btrue(:,3),'b');hold on; grid on;
plot(t_mag,Best(:,3),'r');
xlabel('Time (sec)')

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
title('Noisy $B$ (mG)','fontsize',font_size,'Interpreter','latex')
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

fig=gcf;
figure(fig.Number+1)
subplot(321);
plot(t_mag,Mag_Clean_S(:,1),'r');
grid;ylabel('$B^{S}_{S,1}$','fontsize',font_size,'Interpreter','latex');
title('Error Free $B$ (mG)','fontsize',font_size,'Interpreter','latex')
subplot(322);
plot(t_mag,CalMag(:,1),'b');
grid;
ylabel('$B^{S}_{S,1}$','fontsize',font_size,'Interpreter','latex');
title('Calibrated $B$ (mG)','fontsize',font_size,'Interpreter','latex')
subplot(323);
plot(t_mag,Mag_Clean_S(:,2),'r');grid;
ylabel('$B^{S}_{S,2}$','fontsize',font_size,'Interpreter','latex');
subplot(324);
plot(t_mag,CalMag(:,2),'b');grid;
ylabel('$B^{S}_{S,2}$','fontsize',font_size,'Interpreter','latex');
subplot(325);
plot(t_mag,Mag_Clean_S(:,3),'r');
grid;
ylabel('$B^{S}_{S,3}$','fontsize',font_size,'Interpreter','latex');
xlabel('Time (sec)','fontsize',font_size,'Interpreter','latex');
subplot(326)
plot(t_mag,CalMag(:,3),'b');
grid;ylabel('$B^{S}_{S,3}$','fontsize',font_size,'Interpreter','latex');
xlabel('Time (sec)','fontsize',font_size,'Interpreter','latex');

% 3D plot of True EMF in S frame %
fig=gcf;
figure(fig.Number+1)
plot3(Mag_Clean_S(:,1),Mag_Clean_S(:,2),Mag_Clean_S(:,3),'bo')
xlabel('$$B^{S}_{S,1}$$','Interpreter','Latex','FontSize',font_size);
ylabel('$$B^{S}_{S,2}$$','Interpreter','Latex','FontSize',font_size);
zlabel('$$B^{S}_{S,3}$$','Interpreter','Latex','FontSize',font_size);
grid;
title('Error Free $B$ (mG)','fontsize',font_size,'Interpreter','latex')

% 3D plot of Corrupted EMF in S frame %
fig=gcf;
figure(fig.Number+1)
plot3(uS(:,1),uS(:,2),uS(:,3),'go')
xlabel('$$B^{S}_{S,1,meas}$$','Interpreter','Latex','FontSize',font_size);
ylabel('$$B^{S}_{S,2,meas}$$','Interpreter','Latex','FontSize',font_size);
zlabel('$$B^{S}_{S,3,meas}$$','Interpreter','Latex','FontSize',font_size);
grid;
title('Corrupted $B$ (mG)','fontsize',font_size,'Interpreter','latex')

% 3D plot of Noisy EMF in S frame %
fig=gcf;
figure(fig.Number+1)
plot3(CalMag(:,1),CalMag(:,2),CalMag(:,3),'ro')
xlabel('$$B^{S}_{S,1,cal}$$','Interpreter','Latex','FontSize',font_size);
ylabel('$$B^{S}_{S,2,cal}$$','Interpreter','Latex','FontSize',font_size);
zlabel('$$B^{S}_{S,3,cal}$$','Interpreter','Latex','FontSize',font_size);
grid;
title('Calibrated $B$ (mG)','fontsize',font_size,'Interpreter','latex')

