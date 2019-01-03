function [xdot,lambda] = VecField(x,u,act01,act10)

%--- Vector field ---------------------------------------------------------
%    Computes the time derivative of the 8x1 state x (4 positions, 4 
%    velocities) for a given current state, active constraints, and 
%    applied input u. The function moreover returns the constraint forces
%    lambda (zero when the constraint is inactive). Validity of these 
%    forces is checked and an error is produced when they are negative.
%    Baumgarte constraint stabilization is incorporated.
%
%    Mark Rijnen
%    TU/e, 28-6-2018
%--------------------------------------------------------------------------

% Parameter
wn    = 20;

% Construct components equations of motion
q     = x(1:4);
qd    = x(5:8);

M = MassMat(q);
h = GenTorq(q,qd);
S = [eye(3);zeros(1,3)];

if act01 == 1 && act10 == 1
    [gamma01,J01,r01] = GrdFunc01(q,qd);
    [gamma10,J10,r10] = GrdFunc10(q,qd);
    gamma = [gamma01;gamma10];
    r     = [r01;r10];
    J     = [J01;J10];    
elseif act01 == 1
    [gamma,J,r] = GrdFunc01(q,qd);
elseif act10 == 1
    [gamma,J,r] = GrdFunc10(q,qd);  
elseif act01 == 0 && act10 == 0
    gamma = [];
    J     = [];
    r     = [];
else
    error('Incorrect value for act01 or act10 (as part of the state).')
end

% Construct and solve DAE
A = [M, -J.';J, zeros(length(gamma))];
if ~isempty(gamma)
    B = [-h+S*u; -wn^2*gamma-2*wn*J*qd-r];
else
    B = -h+S*u;
end
xi = A\B;
qdd = xi(1:4);
lambda = xi(5:end);

if isempty(lambda)
    lambda = [0;0];
else
    lambda = [act01==1;act10==1].*lambda;
end

% Check contraint forces
if min(lambda) < 0
    warning('Invalid constraint force encountered.')
    lambda
    act01
    act10
end

% Time derivative of state
xdot = [qd;qdd];

end

