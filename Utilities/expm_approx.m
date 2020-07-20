function [out] = expm_approx(F,dt)

n = length(F);
out = eye(n) + dt*F + 0.5 * (dt*F)^2;

end