function h = GenTorq(q,qd)

%--- Generalized torques --------------------------------------------------
%    Computes the part of the second order dynamics other than inertia,   
%    actuation, and constraint contributions. In particular, the column   
%    h contains Coriolis forces and the spring and damper contributions   
%    of the door (gravity is neglected as a horizontal configuration is   
%    assumed.
%
%    Mark Rijnen
%    TU/e, 28-6-2018
%--------------------------------------------------------------------------

global m2 m3 L1 L2 L3 ck cd 

            qd1 = qd(1);
q2 = q(2);  qd2 = qd(2);
q3 = q(3);  qd3 = qd(3);
q4 = q(4);  qd4 = qd(4);

h1 = -L1*L2*sin(q2)*(m2/2+m3)*qd2^2 - L2*L3*m3*sin(q3)*(qd3^2/2+qd1*qd3+qd2*qd3) - L1*L2*(m2+2*m3)*sin(q2)*qd1*qd2 - 0.5*L1*L3*m3*sin(q2+q3)*(qd2^2+qd3^2+2*qd1*qd3+2*qd1*qd2+2*qd2*qd3);
h2 =  L1*L2*sin(q2)*(m2/2+m3)*qd1^2 - L2*L3*m3*sin(q3)*(qd3^2/2+qd1*qd3+qd2*qd3) + 0.5*L1*L3*m3*sin(q2+q3)*qd1^2;
h3 = 0.5*L3*m3*( L1*sin(q2+q3)*qd1^2 + L2*sin(q3)*(qd1+qd2)^2 );
h4 = cd*qd4 + ck*q4;

h  = [h1;h2;h3;h4];

end