function [qR,qdR,qddR] = InvKin(Y,Yd,Ydd)

%--- Inverse kinematics ---------------------------------------------------
%    Computes the joint angles of the planar RRR robot arm given a        
%    desired position (Y(1),Y(2)) and orientation (Y(3)) of the end-effector.
%    The function also computes the joint velocity and acceleration if the 
%    velocities and accelerations of the end-effector are provided.
%
%    Mark Rijnen
%    TU/e, 28-6-2018
%--------------------------------------------------------------------------

global L1 L2 L3

% Position
x = Y(1);   y = Y(2);   th = Y(3);

o2x = x - L3*cos(th);
o2y = y - L3*sin(th);
D   = (o2x^2 + o2y^2 - L1^2 - L2^2)/(2*L1*L2);
xtest = [];
if abs(D) > 1
    if isempty(xtest)
        xtest = x;
        ytest = y;
        thtest = th;
    D = NaN;
    end
end

q1 = atan2(o2y,o2x)-atan2(L2*sqrt(1-D^2),L1+L2*D);
q2 = atan2(sqrt(1-D^2),D);
q3 = th-atan2(o2y,o2x)+atan2(L2*sqrt(1-D^2),L1+L2*D)-atan2(sqrt(1-D^2),D);

qR = [q1;q2;q3];

% Velocity   
if nargin == 2 && (nargout == 1 || nargout == 2)
    q12 = q1+q2;        q123 = q12+q3;
    J = [-L1*sin(q1)-L2*sin(q12)-L3*sin(q123), -L2*sin(q12)-L3*sin(q123), -L3*sin(q123);
          L1*cos(q1)+L2*cos(q12)+L3*cos(q123),  L2*cos(q12)+L3*cos(q123),  L3*cos(q123);
                                            1,                         1,             1];
    qdR = J\Yd;
    
% Velocity and acceleration    
elseif nargin == 3 
    q12 = q1+q2;        q123 = q12+q3;
    J = [-L1*sin(q1)-L2*sin(q12)-L3*sin(q123), -L2*sin(q12)-L3*sin(q123), -L3*sin(q123);
          L1*cos(q1)+L2*cos(q12)+L3*cos(q123),  L2*cos(q12)+L3*cos(q123),  L3*cos(q123);
                                            1,                         1,             1];
    qdR = J\Yd;
    
    qd1  = qdR(1);      qd12 = qdR(1)+qdR(2);       qd123 = qd12+qdR(3);
    r = [-L1*qd1^2*cos(q1)-L2*qd12^2*cos(q12)-L3*qd123^2*cos(q123);
         -L1*qd1^2*sin(q1)-L2*qd12^2*sin(q12)-L3*qd123^2*sin(q123);
                                                                 0];
    qddR = J\(Ydd-r);
    
elseif nargin ~= 1 || nargout ~= 1
    error('Incorrect number of inputs or outputs to InvKin.m.')
end


end