syms q0 q1 q2 q3
assume(q0,'real')
assume(q1,'real')
assume(q2,'real')
assume(q3,'real')
phi=atan2(2*(q0*q1+q2*q3),1-2*(q1^2+q2^2));
the=asin(2*(q0*q2-q3*q1));
psi=atan2(2*(q0*q3+q1*q2),1-2*(q2^2+q3^2));
dpsi_dq0 = diff(psi,q0);
dpsi_dq1 = diff(psi,q1);
dpsi_dq2 = diff(psi,q2);
dpsi_dq3 = diff(psi,q3);
dtheta_dq0 = diff(the,q0);
dtheta_dq1 = diff(the,q1);
dtheta_dq2 = diff(the,q2);
dtheta_dq3 = diff(the,q3);
dphi_dq0 = diff(phi,q0);
dphi_dq1 = diff(phi,q1);
dphi_dq2 = diff(phi,q2);
dphi_dq3 = diff(phi,q3);

G = [dpsi_dq0 dpsi_dq1 dpsi_dq2 dpsi_dq3;
     dtheta_dq0 dtheta_dq1 dtheta_dq2 dtheta_dq3;
     dphi_dq0 dphi_dq1 dphi_dq2 dphi_dq3];
G(1,1)
G(1,2)
G(1,3)
G(1,4)
G(2,1)
G(2,2)
G(2,3)
G(2,4)
G(3,1)
G(3,2)
G(3,3)
G(3,4)