function [gamma10,J10,r10] = GrdFunc10(q,qd)

%--- Guard function 2 -----------------------------------------------------
%    Computes the contact distance of the right end-effector contact      
%    point that we indicate as contact 2 (10 in binary). It moreover      
%    gives the related Jacobian as well as the term r=Jdot*qdot that      
%    appears in the second time derivative of the guard function.      
%
%    Mark Rijnen
%    TU/e, 28-6-2018
%--------------------------------------------------------------------------

global L1 L2 L3 w3 w4 d3 Dx Dy 
% syms L1 L2 L3 w3 w4 d3 Dx Dy 
% syms q1 q2 q3 q4 qd1 qd2 qd3 qd4
% 
% qd = [qd1;qd2;qd3;qd4];
% q12  = q1+q2;
% q123 = q1+q2+q3;
q1   = q(1);
q12  = q(1)+q(2);
q123 = q(1)+q(2)+q(3);
q4   = q(4);

gamma10 = Dx*sin(q4)+Dy*cos(q4)-L1*sin(q1-q4)-L2*sin(q12-q4)-L3*sin(q123-q4)+(w3-d3)*cos(q123-q4)-w4/2;

a = L1*cos(q1-q4);
b = L2*cos(q12-q4);
c = L3*cos(q123-q4);
d = -d3*sin(q123-q4);
e = -w3*sin(q123-q4);

da = -L1*sin(q1-q4);
db = -L2*sin(q12-q4);
dc = -L3*sin(q123-q4);
dd = (w3-d3)*cos(q123-q4);

J10 = [-a-b-c-d+e, -b-c-d+e, -c-d+e, a+b+c+d-e+Dx*cos(q4)-Dy*sin(q4)];
R10 = [-da-db-dc-dd, -db-dc-dd, -dc-dd,  da+db+dc+dd;
          -db-dc-dd, -db-dc-dd, -dc-dd,     db+dc+dd;
             -dc-dd,    -dc-dd, -dc-dd,        dc+dd;
        da+db+dc+dd,  db+dc+dd,  dc+dd, -da-db-dc-dd-Dx*sin(q4)-Dy*cos(q4)];      
             
r10 = qd.'*R10*qd;
% r10test = L1*(qd1-qd4)^2*sin(q1-q4) + L2*(qd1+qd2-qd4)^2*sin(q12-q4) + L3*(qd1+qd2+qd3-qd4)^2*sin(q123-q4) ...
%     + d3*(qd1+qd2+qd3-qd4)^2*cos(q123-q4) - Dx*qd4^2*sin(q4) - Dy*qd4^2*cos(q4) -w3*(qd1+qd2+qd3-qd4)^2*cos(q123-q4);
end