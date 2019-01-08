clear; close all; clc;

% Parameters
load ref_motion.mat;
param;
t0  = 0;                                % initial time
dt  = 0.01;                             % time step
T   = t_ref{end}(end);                  % final time
x0  = [q_ref{1}(1,:).';qd_ref{1}(1,:).']; % initial condition for reference
Tol = 1e-6;

K   = {8*[1^2*eye(3),[0;0;0],2*1*eye(3),[0;0;0]], ...
       [10^2*[1 0 0 0;0 1 0 0;0 0 0 0],2*10*[1 0 0 0;0 1 0 0;0 0 0 0]]};
erange = 0.8 + rand(1);
x0e = x0 + 1*[-0.05;-0.03;-0.02;0.00; -0.02;-0.05;0.03;0.00];

%**************** Define the hybrid system *******************************%

hybsys.Nm = 4;                          % number of modes
hybsys.Ns = [8,8,8,8];                  % number of states per mode
hybsys.Ni = [3,3,3,3];                  % number of inputs per mode

f1 = @(x,u,~,~) VecField(x,u,0,0);      % vector field mode 1
f2 = @(x,u,~,~) VecField(x,u,1,0);      % vector field mode 2
f3 = @(x,u,~,~) VecField(x,u,0,1);      % vector field mode 3
f4 = @(x,u,~,~) VecField(x,u,1,1);      % vector field mode 4
hybsys.f  = {f1,f2,f3,f4};              % define hybrid system vector fields

c1 = @(x,~,~,~) min(GrdFunc01(x(1:4),x(5:8)),GrdFunc10(x(1:4),x(5:8))); % guard for mode 1
c2 = @(x,~,~,~) GrdFunc10(x(1:4),x(5:8));                               % guard for mode 2
c3 = @(x,~,~,~) GrdFunc01(x(1:4),x(5:8));                               % guard for mode 3
c4 = @(x,~,~,~) 1;                                                      % guard for mode 4
hybsys.c  = {c1,c2,c3,c4};              % define hybrid system guard conditions

e12 = @(x,~,~,~) GrdFunc01(x(1:4),x(5:8)) <= Tol && GrdFunc10(x(1:4),x(5:8)) > Tol;
e13 = @(x,~,~,~) GrdFunc01(x(1:4),x(5:8)) > Tol  && GrdFunc10(x(1:4),x(5:8)) <= Tol;
e14 = @(x,~,~,~) GrdFunc01(x(1:4),x(5:8)) <= Tol && GrdFunc10(x(1:4),x(5:8)) <= Tol;
e24 = @(x,~,~,~) GrdFunc10(x(1:4),x(5:8)) <= Tol;
e34 = @(x,~,~,~) GrdFunc01(x(1:4),x(5:8)) <= Tol;
hybsys.e  = {[], e12, e13, e14; 
             [], [],  [],  e24; 
             [], [],  [],  e34; 
             [], [],  [],  []};         % define hybrid system boolean edge functions

g12 = @(x,~,~) ImpMap(x,1,0);           % reset for transition mode 1 - mode 2
g13 = @(x,~,~) ImpMap(x,0,1);           % reset for transition mode 1 - mode 3
g14 = @(x,~,~) ImpMap(x,1,1);           % reset for transition mode 1 - mode 4
g24 = @(x,~,~) ImpMap(x,1,1);           % reset for transition mode 2 - mode 4
g34 = @(x,~,~) ImpMap(x,1,1);           % reset for transition mode 3 - mode 4                             % reset for transition mode 2 - mode 1
hybsys.g  = {[], g12, g13, g14; 
             [], [],  [],  g24; 
             [], [],  [],  g34; 
             [], [],  [],  []};         % define hybrid system reset maps

%*************************************************************************%

% Define the input (for creating the reference)
Ka = 0*[1^2*eye(3),[0;0;0],2*1*eye(3),[0;0;0]];
ff = {@(t,x,~) interp1(t_ref{1},mu{1},t).', ...
      @(t,x,~) interp1(t_ref{1},mu{1},t).', ...
      @(t,x,~) interp1(t_ref{1},mu{1},t).', ...
      @(t,x,~) interp1(t_ref{2},mu{2},t).'};         % control input for each mode

cntrl = {@(t,x,~) interp1(t_ref{1},mu{1},t).'-Ka*(x-interp1(t_ref{1},[q_ref{1},qd_ref{1}],t).'), ...
         @(t,x,~) interp1(t_ref{1},mu{1},t).'-Ka*(x-interp1(t_ref{1},[q_ref{1},qd_ref{1}],t).'), ...
         @(t,x,~) interp1(t_ref{1},mu{1},t).'-Ka*(x-interp1(t_ref{1},[q_ref{1},qd_ref{1}],t).'), ...
         @(t,x,~) interp1(t_ref{2},mu{2},t).'-Ka*(x-interp1(t_ref{2},[q_ref{2},qd_ref{2}],t).')};         % control input for each mode
      
% Simulate the system
ref_traj = HybridSim(hybsys,ff,t0:dt:T,x0,1);
optionsExt.dom = 'T';
optionsExt.T   = 0.2;
optionsExt.dt  = dt;
cntrl{4} = @(t,x,~) interp1(t_ref{2},mu{2},t).'+Ka*(x-interp1(t_ref{2},[q_ref{2},qd_ref{2}],t).');
ref_traj_ext = ExtndHybTraj(ref_traj,hybsys,cntrl,optionsExt);

% Simulate tracking of the reference
mu           = @(t,s) interp1(ref_traj_ext.t{s},ref_traj_ext.u{s},t).';                     % feedforward
alpha        = @(t,s) interp1(ref_traj_ext.t{s},ref_traj_ext.x{s}(:,1:end-1),t).';          % reference state
cntrl_trc    = {@(t,x,~) mu(t,1) - K{1}*(x-alpha(t,1)), ...
                @(t,x,~) diag([1,1,1])*mu(t,1), ...
                @(t,x,~) diag([1,1,1])*mu(t,1), ...
                @(t,x,~) mu(t,2) - K{2}*(x-alpha(t,2))};   % controller (using pushing sequence) 
traj_trc     = HybridSim(hybsys,cntrl_trc,t0:dt:T,x0e,1);

tue = tue_color;
options.height          = 14;
options.width           = 11;
options.evntlines       = 1;
options.grid            = 2;
options.cntr            = 'num';
options.barcolor        = [tue.db;tue.g;tue.c;tue.b];
options.marker          = {'.','.';'.','.';'.','.'};
options.markersize      = [7*ones(3,1),7*ones(3,1)];
options.fontsize.axes   = 12;
options.fontsize.labels = 13;
options.fontsize.text   = 13;
options.spacingy        = 0.02;
%% Plot reference traj
options.linecolor       = {tue.r,tue.r};
options.linestyle       = {':','-'};
options.cntrbar         = [0,1];

options.labels          = {'$q_1$ [rad]','$q_2$ [rad]','$q_3$ [rad]'};
options.legend          = {'$\bar{\mathbf{\alpha}}$','$\mathbf{\alpha}$'};
signals.x = [1,2,3];
signals.u = [];
grd = [3,1];
figure
PlotHybTraj({ref_traj_ext,ref_traj},signals,grd,options)
movegui('northwest')

options.labels          = {'$\dot{q}_1$ [rad/s]','$\dot{q}_2$ [rad/s]','$\dot{q}_3$ [rad/s]'};
options.legend          = {'$\bar{\mathbf{\alpha}}$','$\mathbf{\alpha}$'};
signals.x = [5,6,7];
signals.u = [];
grd = [3,1];
figure
PlotHybTraj({ref_traj_ext,ref_traj},signals,grd,options)
movegui('north')
% 
options.labels          = {'$\mu_1$ [Nm]','$\mu_2$ [Nm]','$\mu_3$ [Nm]'};
options.legend          = {'$\bar{\mathbf{\mu}}$','$\mathbf{\mu}$'};
signals.x = [];
signals.u = [1,2,3];
grd = [3,1];
figure
PlotHybTraj({ref_traj_ext,ref_traj},signals,grd,options)
movegui('northeast')

% AnimSys(ref_traj.t,ref_traj.x,1)

%% Plot perturbed traj
options.linecolor       = {tue.r};
options.linestyle       = {'-'};
options.cntrbar         = [1];

options.labels          = {'$q_1$ [rad]','$q_2$ [rad]','$q_3$ [rad]'};
options.legend          = {'$\mathbf{x}^{\epsilon}$'};
signals.x = [1,2,3];
signals.u = [];
grd = [3,1];
figure
PlotHybTraj({traj_trc},signals,grd,options)
movegui('northwest')

options.labels          = {'$\dot{q}_1$ [rad/s]','$\dot{q}_2$ [rad/s]','$\dot{q}_3$ [rad/s]'};
options.legend          = {'$\mathbf{x}^{\epsilon}$'};
signals.x = [5,6,7];
signals.u = [];
grd = [3,1];
figure
PlotHybTraj({traj_trc},signals,grd,options)
movegui('north')

options.labels          = {'$u_1$ [Nm]','$u_2$ [Nm]','$u_3$ [Nm]'};
options.legend          = {'$\mathbf{u}$'};
signals.x = [];
signals.u = [1,2,3];
grd = [3,1];
figure
PlotHybTraj({traj_trc},signals,grd,options)
movegui('northeast')

AnimSys(traj_trc.t,traj_trc.x,1)

%% Plot reference + perturbed traj
% AnimSysWRef(ref_traj.t,ref_traj.x,traj_trc.t,traj_trc.x,1)

% 
% Compute constraint forces
figure
for i = 2:3
    lambda = zeros(2,length(traj_trc.t{i}));
    if traj_trc.m(i) == 1
        for j = 1:length(traj_trc.t{i})
            [~,lambda(:,j)] = VecField(traj_trc.x{i}(j,1:8).',traj_trc.u{i}(j,:).',0,0);
        end
    elseif traj_trc.m(i) == 2
        for j = 1:length(traj_trc.t{i})
            [~,lambda(:,j)] = VecField(traj_trc.x{i}(j,1:8).',traj_trc.u{i}(j,:).',1,0);
        end
    elseif traj_trc.m(i) == 3
        for j = 1:length(traj_trc.t{i})
            [~,lambda(:,j)] = VecField(traj_trc.x{i}(j,1:8).',traj_trc.u{i}(j,:).',0,1);
        end
    else
        for j = 1:length(traj_trc.t{i})
            [~,lambda(:,j)] = VecField(traj_trc.x{i}(j,1:8).',traj_trc.u{i}(j,:).',1,1);
        end
    end
    plot(traj_trc.t{i},lambda(1,:).','Color',tue.r,'LineWidth',2)
    hold on
    plot(traj_trc.t{i},lambda(2,:).','Color',tue.db,'LineWidth',2)
end
grid minor
legend('$\lambda_{n,1}$','$\lambda_{n,2}$','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\lambda_n$ [N]','Interpreter','latex')

figure
for i = 2:2
    lambda = zeros(2,length(ref_traj.t{i}));
    if ref_traj.m(i) == 1
        for j = 1:length(ref_traj.t{i})
            [~,lambda(:,j)] = VecField(ref_traj.x{i}(j,1:8).',ref_traj.u{i}(j,:).',0,0);
        end
    elseif ref_traj.m(i) == 2
        for j = 1:length(ref_traj.t{i})
            [~,lambda(:,j)] = VecField(ref_traj.x{i}(j,1:8).',ref_traj.u{i}(j,:).',1,0);
        end
    elseif ref_traj.m(i) == 3
        for j = 1:length(ref_traj.t{i})
            [~,lambda(:,j)] = VecField(ref_traj.x{i}(j,1:8).',ref_traj.u{i}(j,:).',0,1);
        end
    else
        for j = 1:length(ref_traj.t{i})
            [~,lambda(:,j)] = VecField(ref_traj.x{i}(j,1:8).',ref_traj.u{i}(j,:).',1,1);
        end
    end
    plot(ref_traj.t{i},lambda(1,:).','Color',tue.r,'LineWidth',2)
    hold on
    plot(ref_traj.t{i},lambda(2,:).','Color',tue.db,'LineWidth',2)
end
grid minor
legend('$\lambda_{n,1}$','$\lambda_{n,2}$','Interpreter','latex')
xlabel('$t$ [s]','Interpreter','latex')
ylabel('$\lambda_n$ [N]','Interpreter','latex')