function  [dot_x] = ODEs(t,x,const)
%ODES  function for integration using ode45 in orbit determination sim.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [dot_x] =ODES(t,x,const) returns x_dot = f(x,t) by specifying the 
% differential equations of the system in first-order form.
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
% Updated 11/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call constants %%
mu1 = const.mu1;
J2 = const.J2;
Re = const.Re;
m_b = const.m_b;
I_b = const.I_b;

eye_3 = [0;0;1]; % third column of the identity matrix

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Dynamics %%
% First, extract translation states in a convenient form. 
r_a = x(1:3);
dot_r_a = x(4:6);

% Norm of position vector
r_norm = sqrt(r_a'*r_a);

% Then extract attitude states in convenient form
q=x(7:10);
eta=q(1);
epsilon=q(2:4);
omega_ba=x(11:13);

% Calculate Gamma matrix
Gamma_b = GammaQuaternion(q);

% Gravitational force
fg_a = -mu1/r_norm^3*r_a;

% Additional gravitational force due to J2 perturbation
fg_J2_a = 3*mu1*J2*Re^2/(2*r_norm^5)*((5*r_a(3)^2/r_norm^2-1)*r_a - 2*r_a(3)*eye_3);

% Calculate magnetic torque
Cba = quat2dcm(q');
r_b = Cba*r_a; % Position in b frame
b_a = EarthMagField(r_a,t); % Earth's magnetic field
b_b = Cba*b_a; % Earth's magnetic field in b frame
if m_b==0
    tau_mag_b=0;             % magnetic torque
else
tau_mag_b = sk(m_b)*b_b; % magnetic torque
end

% Gravity gradient torque
tau_gg_b = 3*mu1/r_norm^5*sk(r_b)*I_b*r_b;

% Total torque
tau_b = tau_gg_b + tau_mag_b;
tau_b = 0;

% Total force
f_a = fg_a + fg_J2_a;

% Angular Momentum
h_ba=x(14:16);

% Form dot_x = f(x,u) system.
ddot_r_a = f_a;
dot_q = Gamma_b*omega_ba;
dot_omega_ba = I_b\(-sk(omega_ba)*I_b*omega_ba + tau_b);
dot_h_ba = -sk(omega_ba)*h_ba+tau_b;

dot_x = [dot_r_a;ddot_r_a;dot_q;dot_omega_ba;dot_h_ba];
end


