clear; close all; clc;

% Load parameters
param

% Set motion parameters
t_vec = [0,1,2];                    % Start time, impact time, end time
dt    = 0.001;                      % Time step
Dt    = 0.5;                        % Extension time
y{1}  = [0.65,  0,        0; 
         0.45,  Dy-w4/2,  pi/2;
         0.35,  1.35*Dy,  pi/1.5];  % Position and orientation of the end-effector at beginning, impact, and end of (extended) phase 1
                                    % The door/plank is assumed to be at rest during this phase 
y{2}  = [0.4, pi/5];                % Position of the end-effector along the door and angle q4 at the beginning and end of (extended) phase 2. The begin position for q4 is not used, but is chosen such that the contact forces remain positive.
                                    % The post-impact position is determined from phase 1 for consistency
v{1}  = [0,     0,     0
         -0.15, 0.45,  -1;
         -0.15, 1.2,  0];           % (Trans. and rot.) velocity of the end-effector at beginning, impact, and end of (extended) phase 1
                                    % The door/plank is assumed to be at rest during this phase
v{2}  = [0,   0];                   % Velocity of the end-effector along the door and dq4/dt at the beginning and end of (extended) phase 2. The begin velocity for q4 is not used, but is chosen such that the contact forces remain positive.
                                    % The post-impact speed is determined from phase 1 using the impact map (for consistency) 

a{1}  = [0,     0,    0;
         -0.03, 2.5,  -.25;
         0,     0,    0];           % (Trans. and rot.) acceleration of the end-effector at beginning, impact, and end of (extended) phase 1
% a{1}(2,2) = 3;
a{2}  = [-2, 0];                    % Acceleration of the end-effector along the door and d2q4/dt2 at the beginning and end of (extended) phase 2. The begin acc for q4 is not used, but is chosen such that the contact forces remain positive.                                    % Continuity of the x component of the 2d acceleration is used to compute the post impact acc. 

qdd4_pi  = 2;                       % Post-impact angular acceleration of the door 
qddd4_pi = 0;                       % Post-impact angular jerk of the door

% Create motion **********************************************************%

% Polynomial shape functions (in order y0,v0,a0,yf,vf,af)
C = @(t,t0,tf) [((t - tf).^3.*(6*t.^2 - 15*t*t0 + 3*t*tf + 10*t0^2 - 5*t0*tf + tf^2))/(t0 - tf)^5, ...
                                          -((t - t0).*(t - tf).^3.*(3*t - 4*t0 + tf))/(t0 - tf)^4, ...
                                                       ((t - t0).^2.*(t - tf).^3)/(2*(t0 - tf)^3), ...
               -((t - t0).^3.*(6*t.^2 + 3*t*t0 - 15*t*tf + t0^2 - 5*t0*tf + 10*tf^2))/(t0 - tf)^5, ...
                                          -((t - t0).^3.*(t - tf).*(3*t + t0 - 4*tf))/(t0 - tf)^4, ...
                                                      -((t - t0).^3.*(t - tf).^2)/(2*(t0 - tf)^3)];

Cd = @(t,t0,tf) [                                           (30*(t - t0).^2.*(t - tf).^2)/(t0 - tf)^5, ...
                 ((t - tf).^2.*(- 15*t.^2 + 28*t*t0 + 2*t*tf - 12*t0^2 - 4*t0*tf + tf^2))/(t0 - tf)^4, ...
                                         -((t - t0).*(t - tf).^2.*(3*t0 - 5*t + 2*tf))/(2*(t0 - tf)^3), ...
                                                           -(30*(t - t0).^2.*(t - tf).^2)/(t0 - tf)^5, ...
                 ((t - t0).^2.*(- 15*t.^2 + 2*t*t0 + 28*t*tf + t0^2 - 4*t0*tf - 12*tf^2))/(t0 - tf)^4, ...
                                         ((t - t0).^2.*(t - tf).*(2*t0 - 5*t + 3*tf))/(2*(t0 - tf)^3)];

Cdd = @(t,t0,tf) [                          -(60*(t - t0).*(t - tf).*(t0 - 2*t + tf))/(t0 - tf)^5, ...
                                         (12*(t - t0).*(t - tf).*(2*t0 - 5*t + 3*tf))/(t0 - tf)^4, ...
                   ((t - tf).*(10*t.^2 - 12*t*t0 - 8*t*tf + 3*t0^2 + 6*t0*tf + tf^2))/(t0 - tf)^3, ...
                                             (60*(t - t0).*(t - tf).*(t0 - 2*t + tf))/(t0 - tf)^5, ...
                                         (12*(t - t0).*(t - tf).*(3*t0 - 5*t + 2*tf))/(t0 - tf)^4, ...
                  -((t - t0).*(10*t.^2 - 8*t*t0 - 12*t*tf + t0^2 + 6*t0*tf + 3*tf^2))/(t0 - tf)^3];
               
% First part of the motion
t{1} = [t_vec(1):dt:t_vec(2),t_vec(2)+dt:dt:t_vec(2)+Dt];
Y    = C((t_vec(1):dt:t_vec(2)).',t_vec(1),t_vec(2))*[y{1}(1,:);v{1}(1,:);a{1}(1,:);y{1}(2,:);v{1}(2,:);a{1}(2,:)];
Yd   = Cd((t_vec(1):dt:t_vec(2)).',t_vec(1),t_vec(2))*[y{1}(1,:);v{1}(1,:);a{1}(1,:);y{1}(2,:);v{1}(2,:);a{1}(2,:)];
Ydd  = Cdd((t_vec(1):dt:t_vec(2)).',t_vec(1),t_vec(2))*[y{1}(1,:);v{1}(1,:);a{1}(1,:);y{1}(2,:);v{1}(2,:);a{1}(2,:)];
Y    = [Y;C((t_vec(2)+dt:dt:t_vec(2)+Dt).',t_vec(2),t_vec(2)+Dt)*[y{1}(2,:);v{1}(2,:);a{1}(2,:);y{1}(3,:);v{1}(3,:);a{1}(3,:)]];
Yd   = [Yd;Cd((t_vec(2)+dt:dt:t_vec(2)+Dt).',t_vec(2),t_vec(2)+Dt)*[y{1}(2,:);v{1}(2,:);a{1}(2,:);y{1}(3,:);v{1}(3,:);a{1}(3,:)]];
Ydd  = [Ydd;Cdd((t_vec(2)+dt:dt:t_vec(2)+Dt).',t_vec(2),t_vec(2)+Dt)*[y{1}(2,:);v{1}(2,:);a{1}(2,:);y{1}(3,:);v{1}(3,:);a{1}(3,:)]];

q{1}  = zeros(length(t{1}),4);      qd{1} = q{1};       qdd{1} = q{1};
for i = 1:length(t{1})
    [q{1}(i,1:3),qd{1}(i,1:3),qdd{1}(i,1:3)] = InvKin(Y(i,:).',Yd(i,:).',Ydd(i,:).');
end

% Post-impact boundary conditions
q_min   = interp1(t{1},q{1},t_vec(2),'pchip').';
qd_min  = interp1(t{1},qd{1},t_vec(2),'pchip').';
qdd_min = interp1(t{1},qdd{1},t_vec(2),'pchip').';
x_plus  = ImpMap([q_min;qd_min],1,1);
qd_plus = x_plus(5:8);
q1  = q_min(1);      q2  = q_min(2);      q3  = q_min(3);      q4  = q_min(4);
qd1 = qd_plus(1);    qd2 = qd_plus(2);    qd3 = qd_plus(3);    qd4 = qd_plus(4);
y_pi = [Dx*cos(q4)-Dy*sin(q4)+L1*cos(q1-q4)+L2*cos(q1+q2-q4)+L3*cos(q1+q2+q3-q4), q4]; %Determine q_D from post-event state
Je4  = [-L1*sin(q1-q4)-L2*sin(q1+q2-q4)-L3*sin(q1+q2+q3-q4), -L2*sin(q1+q2-q4)-L3*sin(q1+q2+q3-q4), -L3*sin(q1+q2+q3-q4), -Dx*sin(q4)-Dy*cos(q4)+L1*sin(q1-q4)+L2*sin(q1+q2-q4)+L3*sin(q1+q2+q3-q4); 0 0 0 1]; %J_D
v_pi = (Je4*qd_plus).';
a_pi = (Je4*qdd_min + [-L1*(qd1-qd4)^2*cos(q1-q4)-L2*(qd1+qd2-qd4)^2*cos(q1+q2-q4)-L3*(qd1+qd2+qd3-qd4)^2*cos(q1+q2+q3-q4)-qd4^2*(Dx*cos(q4)-Dy*sin(q4));0]).';
a_pi(2) = qdd4_pi;

% Second part of the motion
t{2} = [t_vec(2):dt:t_vec(3)];
% Y2    = C((t_vec(2)-Dt:dt:t_vec(2)-dt).',t_vec(2)-Dt,t_vec(2))*[y{2}(1,:);v{2}(1,:);a{2}(1,:);y_pi;v_pi;a_pi];
% Yd2   = Cd((t_vec(2)-Dt:dt:t_vec(2)-dt).',t_vec(2)-Dt,t_vec(2))*[y{2}(1,:);v{2}(1,:);a{2}(1,:);y_pi;v_pi;a_pi];
% Ydd2  = Cdd((t_vec(2)-Dt:dt:t_vec(2)-dt).',t_vec(2)-Dt,t_vec(2))*[y{2}(1,:);v{2}(1,:);a{2}(1,:);y_pi;v_pi;a_pi];
% t_part = (-Dt:dt:-dt).';
% Y2(:,2)   = [ones(length(t_part),1), t_part,0.5*t_part.^2,1/6*t_part.^3]*[y_pi(2);v_pi(2);qdd4_pi;qddd4_pi]; 
% Yd2(:,2)  = [zeros(length(t_part),1),ones(length(t_part),1),t_part,1/2*t_part.^2]*[y_pi(2);v_pi(2);qdd4_pi;qddd4_pi];
% Ydd2(:,2) = [zeros(length(t_part),1),zeros(length(t_part),1),ones(length(t_part),1),t_part]*[y_pi(2);v_pi(2);qdd4_pi;qddd4_pi];

Y2    = [C((t_vec(2):dt:t_vec(3)).',t_vec(2),t_vec(3))*[y_pi;v_pi;a_pi;y{2}(1,:);v{2}(1,:);a{2}(1,:)]]; 
Yd2   = [Cd((t_vec(2):dt:t_vec(3)).',t_vec(2),t_vec(3))*[y_pi;v_pi;a_pi;y{2}(1,:);v{2}(1,:);a{2}(1,:)]];
Ydd2  = [Cdd((t_vec(2):dt:t_vec(3)).',t_vec(2),t_vec(3))*[y_pi;v_pi;a_pi;y{2}(1,:);v{2}(1,:);a{2}(1,:)]];

q{2}  = zeros(length(t{2}),4);      qd{2} = q{2};       qdd{2} = q{2};
for i = 1:length(t{2})
    q4 = Y2(i,2); %obtain q4 from qD
    q{2}(i,:) = [InvKin([cos(q4)*Y2(i,1)+w4/2*sin(q4)-Dx; sin(q4)*Y2(i,1)-w4/2*cos(q4)+Dy; q4+pi/2]).',q4]; %obtain @
    [~,J01,~] = GrdFunc01(q{2}(i,:).',NaN*ones(4,1));
    [~,J10,~] = GrdFunc10(q{2}(i,:).',NaN*ones(4,1));
    q1 = q{2}(i,1);     q2 = q{2}(i,2);     q3 = q{2}(i,3);
    Je4  = [-L1*sin(q1-q4)-L2*sin(q1+q2-q4)-L3*sin(q1+q2+q3-q4), -L2*sin(q1+q2-q4)-L3*sin(q1+q2+q3-q4), -L3*sin(q1+q2+q3-q4), -Dx*sin(q4)-Dy*cos(q4)+L1*sin(q1-q4)+L2*sin(q1+q2-q4)+L3*sin(q1+q2+q3-q4); 0 0 0 1];
    J    = [Je4;J01;J10];
    qd{2}(i,:) = (J\[Yd2(i,:).';0;0]).';  
    [~,~,r01] = GrdFunc01(q{2}(i,:).',qd{2}(i,:).');
    [~,~,r10] = GrdFunc10(q{2}(i,:).',qd{2}(i,:).');
    qd1 = qd{2}(i,1);     qd2 = qd{2}(i,2);     qd3 = qd{2}(i,3);       qd4 = qd{2}(i,4);
    re4 = [-L1*(qd1-qd4)^2*cos(q1-q4)-L2*(qd1+qd2-qd4)^2*cos(q1+q2-q4)-L3*(qd1+qd2+qd3-qd4)^2*cos(q1+q2+q3-q4)-qd4^2*(Dx*cos(q4)-Dy*sin(q4));0];
    qdd{2}(i,:) = (J\([-re4+Ydd2(i,:).';-r01;-r10])).';
end

% Input torques and contact forces
mu{1}     = zeros(length(t{1}),3);      mu{2}     = zeros(length(t{2}),3);
lambda{1} = zeros(length(t{1}),2);      lambda{2} = zeros(length(t{2}),2); 
for i = 1:length(t{1})
    [mu_i,lambda_i] = InvDyn(q{1}(i,:).',qd{1}(i,:).',qdd{1}(i,:).',0,0);
    mu{1}(i,:)      = mu_i.';
    lambda{1}(i,:)  = lambda_i.';
end
for i = 1:length(t{2})
    [mu_i,lambda_i] = InvDyn(q{2}(i,:).',qd{2}(i,:).',qdd{2}(i,:).',1,1);
    mu{2}(i,:)      = mu_i.';
    lambda{2}(i,:)  = lambda_i.';
end

% Check if flow would be possible if one point impacts
x_plus   = ImpMap([q_min;qd_min],1,0);
[~,test01] = VecField(x_plus,mu{1}(t{1}==t_vec(2),:).',1,0)
x_plus   = ImpMap([q_min;qd_min],0,1);
[~,test10] = VecField(x_plus,mu{1}(t{1}==t_vec(2),:).',0,1)

% Save data
t_ref   = t;    q_ref   = q;    qd_ref  = qd;   qdd_ref = qdd;
save ref_motion.mat t_ref q_ref qd_ref qdd_ref mu lambda
clearvars t_ref q_ref qd_ref qdd_ref

% Plot results
figure
plot(t{1},Y(:,1),t{1},Yd(:,1),t{1},Ydd(:,1))
hold on
plot([1 1],[-1 1],'--k')
title('$x$','Interpreter','Latex')
xlabel('$t$','Interpreter','Latex')
legend('pos','vel','acc')

figure
plot(t{1},Y(:,2),t{1},Yd(:,2),t{1},Ydd(:,2))
hold on
plot([1 1],[-6 10],'--k')
title('$y$','Interpreter','Latex')
xlabel('$t$','Interpreter','Latex')
legend('pos','vel','acc')

figure
plot(t{1},Y(:,3),t{1},Yd(:,3),t{1},Ydd(:,3))
hold on
plot([1 1],[-20 20],'--k')
title('$\theta$','Interpreter','Latex')
xlabel('$t$','Interpreter','Latex')
legend('pos','vel','acc')

figure
plot(t{2},Y2(:,1),t{2},Yd2(:,1),t{2},Ydd2(:,1))
title('$o_{3,x}^{4}$','Interpreter','Latex')
xlabel('$t$','Interpreter','Latex')
legend('pos','vel','acc')

figure
plot(t{2},Y2(:,2),t{2},Yd2(:,2),t{2},Ydd2(:,2))
title('$q_4$','Interpreter','Latex')
xlabel('$t$','Interpreter','Latex')
legend('pos','vel','acc')

% figure
% ax(1) = subplot(3,1,1);
% plot(t{1},mu{1}(:,1))
% % hold on
% % plot(t{2},mu{2}(:,1))
% ylabel('$\mu_1$ [Nm]','Interpreter','Latex')
% 
% ax(2) = subplot(3,1,2);
% plot(t{1},mu{1}(:,2))
% % hold on
% % plot(t{2},mu{2}(:,2))
% ylabel('$\mu_2$ [Nm]','Interpreter','Latex')
% 
% ax(3) = subplot(3,1,3);
% plot(t{1},mu{1}(:,3))
% % hold on
% % plot(t{2},mu{2}(:,3))
% ylabel('$\mu_3$ [Nm]','Interpreter','Latex')
% xlabel('$t$ [s]','Interpreter','Latex')
% linkaxes(ax,'x')
% 
% figure
% ax2(1) = subplot(2,1,1);
% plot(t{1},lambda{1}(:,1))
% % hold on
% % plot(t{2},lambda{2}(:,1))
% ylabel('$\lambda_1$ [N]','Interpreter','Latex')
% 
% ax2(2) = subplot(2,1,2);
% plot(t{1},lambda{1}(:,2))
% % hold on
% % plot(t{2},lambda{2}(:,2))
% ylabel('$\lambda_2$ [N]','Interpreter','Latex')
% xlabel('$t$ [s]','Interpreter','Latex')
% linkaxes(ax2,'x')
