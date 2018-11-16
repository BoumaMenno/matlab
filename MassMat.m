function M = MassMat(q)

%--- Mass matrix ----------------------------------------------------------
%    Computes the 4x4 mass matrix for the planar RRR robot (with door) as 
%    a function of the 4x1 generalized coordinates column q 
%
%    Mark Rijnen
%    TU/e, 28-6-2018
%--------------------------------------------------------------------------

global m2 m3 I1o I2o I3o I4o L1 L2 L3 c3x c3y

q2 = q(2);
q3 = q(3);

M11 = I1o + I2o + I3o + L1*L2*(m2+2*m3)*cos(q2) + L2*(L3-c3x)*m3*cos(q3) + 2*L1*(L3-c3x)*m3*cos(q2+q3) + 2*L1*c3y*m3*sin(q2+q3) + 2*L2*c3y*m3*sin(q3);
M12 = I2o + I3o + L1*L2*(m2/2+m3)*cos(q2) + 2*L2*(L3-c3x)*m3*cos(q3) + L1*(L3-c3x)*m3*cos(q2+q3) + L1*c3y*m3*sin(q2+q3) + L2*c3y*m3*sin(q3);
M13 = I3o + L2*(L3-c3x)*m3*cos(q3) + L1*(L3-c3x)*m3*cos(q2+q3) + L1*c3y*m3*sin(q2+q3) + L2*c3y*m3*sin(q3);
M22 = I2o + I3o + 2*L2*(L3-c3x)*m3*cos(q3) + 2*L2*c3y*m3*sin(q3);
M23 = I3o + L2*(L3-c3x)*m3*cos(q3) + L2*c3y*m3*sin(q3);
M33 = I3o;
M44 = I4o;

M = [M11, M12, M13, 0;
     M12, M22, M23, 0;
     M13, M23, M33, 0;
     0,   0,   0,   M44];
 
end