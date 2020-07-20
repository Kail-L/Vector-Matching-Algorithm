function [Cd]=disrw(F,G,T,Rwpsd)
%DISRW  Computes discrete equivalent of continuous noise.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [CD] = disrw(F,G,T,RWPSD) This function discretizes continous noise given
% system dynamics matrices in addition to the state process noise matrix.
% It is essentially performing the integral portion of the discretized
% state propogation equation. 
% 
% SOURCES:
% n/a (Taken from Demoz Gebre toolbox)
%
% INPUT PARAMETERS:
% F = m x n state dynamics matrix, where n is the number of states for our system.
% G = n x m state input matrix.
% T = Discretizing time step (delta T between samples)
% Rwpsd = state process noise covariance matrix
%
% OUTPUT PARAMETERS:
% Cd = discretized input matrix
%
% VARIABLES (note, I am not absolutely sure what each variable really is for):
% ZF = m x n zero matrix
% n = row dimension of G
% m = column dimension of G
% ME = updated matrix for matrix exponential equation
% phi = matrix exponential of ME and T
% phi12 = first row, second column values of Phi
% phi22 = second row, second column values of phi
%
% Kail Laughlin, taken from Demoz Gebre
% Updated 11/25/2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ZF=zeros(size(F));
[n,m]=size(G);
ME=[-F G*Rwpsd*G';
ZF F'];
phi=expm(ME*T);
phi12=phi(1:n,n+1:2*n);
phi22=phi(n+1:2*n,n+1:2*n);
Cd=phi22'*phi12;

