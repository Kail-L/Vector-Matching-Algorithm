function qnorm = quatnorm(q)
q0=q(1);
q1=q(2);
q2=q(3);
q3=q(4);
qnorm=sqrt(q0^2+q1^2+q2^2+q3^2);
end