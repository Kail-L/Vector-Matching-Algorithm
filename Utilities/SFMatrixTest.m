syms q0 q1 q2 q3 wm1 wm2 wm3 SF1 SF2 SF3 dSF1 dSF2 dSF3
assume(q0, 'real')
assume(q1, 'real')
assume(q2, 'real')
assume(q3, 'real')
assume(wm1, 'real')
assume(wm2, 'real')
assume(wm3, 'real')

% assume
SF=[SF1 0 0;
     0 SF2 0;
     0 0 SF3];
 SF1h=SF1+dSF1;
 SF2h=SF2+dSF2;
 SF3h=SF3+dSF3;
 
SFh=[1+SF1h 0 0;
     0 1+SF2h 0;
     0 0 1+SF3h];
iSFh = inv(SFh)
iSFh*SF
q = [q0 q1 q2 q3]';
qv = q(2:4);
w = [wm1; wm2; wm3];
O=1/2*[0 -w';w -sk(w)]*q
iO = quatInverse(O')
QT = [0 -w'; w -sk(w)]*q
QP = quatmult([0;w],q)
eul = [60; 22; 10]*pi/180;
C = eul2dcm(eul);
w = [0.3;0.1;0.5];
D = diag([1.1;1.2;1.3])*C*w;
V = diag(C*w)*[1.1;1.2;1.3];
