function [u,lambda] = InvDyn(q,qd,qdd,act01,act10)

%--- Inverse dynamics -----------------------------------------------------
%    Computes the input u and constraint forces lambda that together
%    achieve a system acceleration qdd for a given position q and velocity
%    qd. The input variables act01 and act10 indicate whether (1) or not 
%    (0) the left, resp. the right constraint are active. Validity of the 
%    constraint forces is checked and an error is produced when one or more 
%    is negative. The pseudoinverse is used in solving the problem when 
%    the equations are under- or overdetermined. In the former, consistency
%    is checked, whereas in the latter situation, the norm of the input and
%    constraint forces is minimized.
%
%    Mark Rijnen
%    TU/e, 28-6-2018
%--------------------------------------------------------------------------

% Construct components equations of motion
M = MassMat(q);
h = GenTorq(q,qd);
S = [eye(3);zeros(1,3)];

if act01 == 1 && act10 == 1
    [~,J01,~] = GrdFunc01(q,qd);
    [~,J10,~] = GrdFunc10(q,qd);
    J     = [J01;J10];    
elseif act01 == 1
    [~,J,~] = GrdFunc01(q,qd);
elseif act10 == 1
    [~,J,~] = GrdFunc10(q,qd);  
elseif act01 == 0 && act10 == 0
    J     = [];
else
    error('Incorrect value for act01 or act10.')
end

% Compute inverse dynamics
if ~isempty(J)
    A  = [S, [sum(J(:,1)),sum(J(:,2)),sum(J(:,3)),sum(J(:,4))].'];
else
    A = S;
end
B  = [M*qdd+h];
xi = A\B;

u      = xi(1:3);
lambda = xi(4:end);

if isempty(lambda)
    lambda = [0;0];
else
    lambda = [act01==1;act10==1]*lambda;
end

% Check consistency if underdetermined system of equations
if size(A,1) > size(A,2)
    if isempty(J)
        err = M*qdd+h-S*u;
    else
        err = M*qdd+h-S*u+J*lambda;
    end
    if norm(err) > 1e-10
        warning('Inconsistency in inverse dynamics detected.')
    end
end

% Check contraint forces
if min(lambda) < 0
    lambda
    warning('Invalid constraint force encountered.')
end

end