function [Cba] = TRIAD(s1_a,s2_a,s1_b,s2_b)
%TRIAD  Calculate Cba using TRIAD algorithm with 2 vector measurements.
% [Cba] = TRIAD(s1_a,s2_a,s1_b,s2_b) solves for Cba using TRIAD algorithm
% using 2 vectors resolved in a and b frames.
%
% INPUT PARAMETERS:
% s1_a = vector s1 resolved in a frame
% s2_a = vector s2 resolved in a frame
% s1_b = vector s1 resolved in b frame
% s2_b = vector s2 resolved in b frame
%
% OUTPUT PARAMETERS:
% Cba = Estimated DCM
%
% Ryan Caverly
% Updated April 2019
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Normalize measurements
s1_hat_a = s1_a/norm(s1_a);
s2_hat_a = s2_a/norm(s2_a);
s1_hat_b = s1_b/norm(s1_b);
s2_hat_b = s2_b/norm(s2_b);

w1_a = s1_hat_a;
w1_b = s1_hat_b;

w2_a = (sk(s1_hat_a)*s2_hat_a)/norm((sk(s1_hat_a)*s2_hat_a));
w2_b = (sk(s1_hat_b)*s2_hat_b)/norm((sk(s1_hat_b)*s2_hat_b));

w3_a = sk(w1_a)*w2_a;
w3_b = sk(w1_b)*w2_b;

Caw = [w1_a w2_a w3_a];
Cbw = [w1_b w2_b w3_b];

Cba = Cbw*Caw';