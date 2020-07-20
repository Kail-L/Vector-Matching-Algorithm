%% Plotting script %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script to plot simulated truth data for a satellite ran through the
% SatAttEst.m code. 

% Kail Laughlin
% 11/25/2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plots.

% Plot of orbit about Earth
EarthPlot(x_out(:,1),x_out(:,2),x_out(:,3),const.Re)


% Plots of positions and velocities
fig=gcf;
figure(fig.Number+1)
subplot(2,1,1) % Plot of x vs time
plot(t,x_out(:,1),'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$x$ (m)','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
subplot(2,1,2) % Plot of x_dot vs time
plot(t,x_out(:,4),'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\dot{x}$ (m/s)','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
saveas(fig,'XXdotVt.png')

fig=gcf;
figure(fig.Number+1)
subplot(2,1,1) % Plot of y vs time
plot(t,x_out(:,2),'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$y$ (m)','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
subplot(2,1,2) % Plot of y_dot vs time
plot(t,x_out(:,5),'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\dot{y}$ (m/s)','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
saveas(fig,'YYdotVt.png')

fig=gcf;
figure(fig.Number+1)
subplot(2,1,1) % Plot of z vs time
plot(t,x_out(:,3),'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$z$ (m)','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
subplot(2,1,2) % Plot of z_dot vs time
plot(t,x_out(:,6),'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\dot{z}$ (m/s)','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
saveas(fig,'ZZdotVt.png')

% Quaternion versus time
fig=gcf;
figure(fig.Number+1)
subplot(4,1,1)
plot(t,q_ba(:,1),'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$q_0$','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
subplot(4,1,2)
plot(t,q_ba(:,2),'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$q_1$','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
subplot(4,1,3)
plot(t,q_ba(:,3),'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$q_2$','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
subplot(4,1,4)
plot(t,q_ba(:,4),'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$q_3$','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)

% 3-2-1 Euler angles versus time
fig=gcf;
figure(fig.Number+1)
subplot(3,1,3)
plot(t,Euler_angles(:,1)*180/pi,'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\phi$ (deg)','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
subplot(3,1,2)
plot(t,Euler_angles(:,2)*180/pi,'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\theta$ (deg)','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
subplot(3,1,1)
plot(t,Euler_angles(:,3)*180/pi,'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\psi$ (deg)','fontsize',font_size,'Interpreter','latex','rotation',0);
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on

% Norm of the quaternion minus 1 (should be zero)
fig=gcf;
figure(fig.Number+1)
plot(t,norm_q,'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$q^Tq - 1$','fontsize',font_size,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on

% Energy, which should be constant
fig=gcf;
figure(fig.Number+1)
subplot(2,1,1) % Plot of total energy vs time (should be constant)
plot(t,E,'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$E$ (J)','fontsize',font_size,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
subplot(2,1,2) % Plot of change in energy vs time (should be zero)
plot(t,E-E(1),'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$E - E(0)$ (J)','fontsize',font_size,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on

% Angular Velocity vs Time
fig=gcf;
figure(fig.Number+1)
subplot(3,1,1)
plot(t,omega_ba(:,1)*180/pi,'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\omega_{ba}^{b1}$ (deg/s)','fontsize',font_size,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
subplot(3,1,2)
plot(t,omega_ba(:,2)*180/pi,'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\omega_{ba}^{b2}$ (deg/s)','fontsize',font_size,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
grid on
subplot(3,1,3)
plot(t,omega_ba(:,3)*180/pi,'Linewidth',line_width);
hold on
xlabel('Time (s)','fontsize',font_size,'Interpreter','latex');
ylabel('$\omega_{ba}^{b3}$ (deg/s)','fontsize',font_size,'Interpreter','latex');
set(gca,'XMinorGrid','off','GridLineStyle','-','FontSize',line_size)
