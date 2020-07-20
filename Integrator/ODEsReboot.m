function  [dot_x] = ODEsReboot(t,x,const)
%ODESREBOOT  RHS of attitude equations used within ODE45 for daily reboot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [dot_x] =ODESREBOOT(t,x,const) returns x_dot = f(x,t) by specifying the 
% differential equations of the system in first-order form. This specific
% ODEs is used to provide a hot-start to the cubesat after a daily reboot
% has occured.
%
% INPUT PARAMETERS:
% t = time
% x = system states
% const = a structure that contains all relevant physical parameters
%
% OUTPUT PARAMETERS:
% dot_x = the first-order differential equation evaluated at x and t
%
% Kail Laughlin
% Updated 3/6/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call constants %%
I_b = const.I_b;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dynamics %%

% Then extract attitude states in convenient form
w=x(1:3);
q=x(4:7);

% Calculate Gamma matrix
Gamma_b = GammaQuaternion(q);

% Form dot_x = f(x,u) system.
dot_omega_ba = I_b\(-sk(w)*I_b*w);
dot_q = Gamma_b*w;

dot_x = [dot_omega_ba;dot_q];
end


