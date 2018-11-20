% Chosen parameters for the planar RRR robot and door (link 4).
% The parameters are (directly) used in AnimSys.m, GenTorq.m, GrdFunc01.m,
% GrdFunc10.m, InvKin.m, and MassMat.m.

global m1 m2 m3 I1o I2o I3o I4o L1 L2 L3 L4 w3 w4 d3 c3x c3y ck cd Dx Dy 

%---- Link 1 -------------------------------------------------------------%
m1 = 4;
L1 = 0.3;
I1 = 0.03;

%---- Link 2 -------------------------------------------------------------%
m2 = 4;
L2 = 0.3;
I2 = 0.03;

%---- Link 3 -------------------------------------------------------------%
m3 = 2;
L3 = 0.1;
I3 = 0.005;
w3 = 0.15;
d3 = 0.03;
c3x = 0.05;
c3y = 0.02;

%---- Link 4 -------------------------------------------------------------%
L4  = 0.8;
w4  = 0.04;
I4o = 5;
ck  = 25;
cd  = 25;

Dx  = 0.1;
Dy  = 0.35;

%---- Effective inertias -------------------------------------------------%
I1o = I1 + 1/4*m1*L1^2 + m2*L1^2 + m3*L1^2;
I2o = I2 + 1/4*m2*L2^2 + m3*L2^2;
I3o = I3 + m3*((L3-c3x)^2 + c3y^2);
