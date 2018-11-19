function x_plus = ImpMap(x_min,act01,act10)

%--- Jump/impact map _-----------------------------------------------------
%    Computes the 8x1 post impact state x_plus given a pre impact state       
%    x_min, knowing which constraints will be active after impact.        
%    The (post impact) activity of the left (resp. right) contact point   
%    is indicated by the variable act01 (resp. act10) where a 1 indicates 
%    that the constraint has been triggered and a 0 means that it is      
%    inactive. Feasibility of the required contact impulses is checked. 
%
%    Mark Rijnen
%    TU/e, 28-6-2018
%--------------------------------------------------------------------------

q_min  = x_min(1:4);
qd_min = x_min(5:8);
M = MassMat(q_min);

if act01 == 1 && act10 == 1
    [~,J01,~] = GrdFunc01(q_min,qd_min);
    [~,J10,~] = GrdFunc10(q_min,qd_min);
    J = [J01;J10];    
elseif act01 == 1
    [~,J01,~] = GrdFunc01(q_min,qd_min);
    J = J01;
elseif act10 == 1
    [~,J10,~] = GrdFunc10(q_min,qd_min);  
    J = J10;
else
    error('Incorrect input to the impact map function ImpMap')
end

xi = [M, -J.'; J, zeros(size(J,1))]\[M*qd_min;zeros(size(J,1),1)];
qd_plus = xi(1:4);
Lambda  = xi(5:end);
if min(Lambda) < 0
    Lambda
    error('Infeasible impact')
end
x_plus = [q_min;qd_plus];

end