function [qhat_ba_reboot,what_ba_reboot,RebootTime] = rebootAtt(qhat_ba,omega_ba,TimeStep,const)
%REBOOTATT function to compute a rough attitude estimate after daily
%reboot has occured on orbit.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [qhat_ba_reboot,what_ba_reboot,RebootTime] =
% REBOOTATT(qhat_ba,omega_ba,const) returns an "updated" attitude state
% after a daily reboot has occured on orbit. To reduce convergence time of
% the EKF, we record the final attitude state (quaternion and bias) before
% a reboot causes data loss. This attitude state is then integrated using
% ODE45 over the timespan of the reboot to provide a hot-start attitude
% condition that is hopefully within 10 degrees of the actual attitude at
% the time of restart. Depending on the TimeStep, the rebooted attitude
% will be more or less accurate (Larger TimeStep = More Accurate). However,
% this will also increase the RebootTime parameter. For the case of this
% simulation, that is not an issue. But when this is implemented on
% hardware, we will need to make this as fast as possible. 
%
% INPUT PARAMETERS:
% qhat_ba = 4x1 or 1x4 quaternion estimat with scalar as first component
% omega_ba = 3x1 or 1x3 body angular rates (rad/s)
% TimeStep = scalar value that signifies the integration dt
% const = a structure that contains all relevant physical parameters
%
% OUTPUT PARAMETERS:
% qhat_ba_reboot = Rebooted 1x4 quaternion 
% omega_ba_reboot = Rebooted 1x3 omega (rad/s)
% RebootTime = How long ODE45 took to run reboot
%
% VARIABLES:
% I_b = moment of inertia matrix taken directly from const structure
% dt = Integration time step
% RL = Length of reboot (s)
% w_prior = Last measured body rates before reboot (rad/s)
% q_prior = Last measured attitude quaternion before reboot
% t_Reboot = Time vector to be input into the ODE45 solver
% IC = Initial condition array to be used in 
%
% Kail Laughlin
% Updated 3/6/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
I_b=const.I_b;
dt=1/TimeStep;
RL=const.RebootLength;
w_prior=[omega_ba(1);
         omega_ba(2);
         omega_ba(3)];
q_prior=[qhat_ba(1);
         qhat_ba(2);
         qhat_ba(3);
         qhat_ba(4)];
     
t_Reboot=1:dt:RL;

IC = [w_prior;q_prior];
options = odeset('AbsTol',1e-9,'RelTol',1e-9); % This changes the integration tolerence.
[t_Reboot,x_Reboot_out] = ode45(@ODEsReboot,t_Reboot,IC,options,const);
qhat_ba_reboot = x_Reboot_out(end,4:7);
what_ba_reboot = x_Reboot_out(end,1:3);
RebootTime=toc;
end