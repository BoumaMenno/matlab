function AnimSys(time,q,showpath,filename)

%--- Animate system -------------------------------------------------------
%    Creates an animation (schematic style) of the 3DOF planar            
%    manipulator interacting with/opening a door. The 4 generalized       
%    coordinates as a function of time are stored in q. Visualization of  
%    the path of the end-effector can be enabled (1) or disabled (0)      
%    using the binary variable 'showpath'. If a string is supplied for    
%    the variable 'filename', an mp4 of the animation is produced.
%
%    Mark Rijnen
%    TU/e, 28-6-2018
%--------------------------------------------------------------------------

global L1 L2 L3 L4 w3 w4 c3x d3 Dx Dy

% Check if a video is requested or not
if nargin == 4
    rec = 1;
elseif nargin == 3
    rec = 0;
else
    error('Incorrect number of input arguments to AnimSys.')
end

% Video parameters and initialization
fps = 30;
dt  = 1/fps;      
if rec == 1
    writerObj = VideoWriter(filename,'MPEG-4');
    writerObj.FrameRate = fps;
    writerObj.Quality   = 100;
    open(writerObj);
end

% Combine data if hybrid trajectory is supplied
if iscell(time) && iscell(q)
    t = [];
    Q = [];
    for i = 1:length(time)
        t = [t;time{i}];
        Q = [Q;q{i}(:,1:4)];
    end
    [time,I,~] = unique(t);
    q = Q(I,:);
    clearvars t Q
elseif iscell(time) || iscell(q)
    error('Incorrect input argument for time or position variable in AnimSys.m')
end

% Transpose q and time if necessary
size1 = size(q,1);  size2 = size(q,2);
if ~(size1 == 4 || size2 == 4)
    error('Dimension of the joint coordinates trajectory vector q is incorrect')
elseif size1 ~= 4
    q = q.';
end
if min(size(time)) ~= 1
    error('Dimension of time vector is incorrect.')
elseif size(time,1) ~= 1
    time = time.';
end

% Resample the trajectory
t    = time(1):dt:time(end);
q1   = interp1(time,q(1,:),t);
q2   = interp1(time,q(2,:),t); 
q3   = interp1(time,q(3,:),t);
q4   = interp1(time,q(4,:),t);
q12  = q1 + q2;
q123 = q1 + q2 + q3;

% Create geometry of the components
r  = 0.04;  lw = 1;  dx = r/2.5;
th = 0:pi/100:pi; 
door  = [0 L4 L4 0 0;  w4/2 w4/2 -w4/2 -w4/2 w4/2; 1 1 1 1 1];
link1 = [0 -L1-r*sin(th) r*sin(th); r r*cos(th) -r*cos(th); ones(1,1+2*length(th))]; 
link2 = [0 -L2-r*sin(th) r*sin(th); r r*cos(th) -r*cos(th); ones(1,1+2*length(th))];
link3 = [[0, 0, -c3x, -L3, -L3];
         [d3, d3-w3, d3-w3, 0, d3];
          ones(1,5)];
circ  = [0.5*r*cos(2*th);0.5*r*sin(2*th)];
hinge = 1.2*[-w4 w4/2*sin(th) -w4 -w4; w4/2 w4/2*cos(th) -w4/2 w4/2];
base  = 1.2*[r r*cos(th) -r r;-2*r r*sin(th) -2*r -2*r];

% Get end-effector position as a function of time
x  = L1*cos(q1)+L2*cos(q12)+L3*cos(q123);
y  = L1*sin(q1)+L2*sin(q12)+L3*sin(q123);

% Color info
cmap = colormap(lines);
col = [1,1,1];

% Animate
% close all;  
figure
movegui('east')
for ii = 1:length(t)
    clf('reset')
    set(gcf,'Color','w')
    
    % Transformation matrices
    H01 = [cos(q1(ii)), -sin(q1(ii)), L1*cos(q1(ii));
           sin(q1(ii)),  cos(q1(ii)), L1*sin(q1(ii));
                     0,            0,              1];
    H02 = [cos(q12(ii)), -sin(q12(ii)), L1*cos(q1(ii))+L2*cos(q12(ii));
           sin(q12(ii)),  cos(q12(ii)), L1*sin(q1(ii))+L2*sin(q12(ii));
                      0,             0,                              1];
    H03 = [cos(q123(ii)), -sin(q123(ii)), L1*cos(q1(ii))+L2*cos(q12(ii))+L3*cos(q123(ii));
           sin(q123(ii)),  cos(q123(ii)), L1*sin(q1(ii))+L2*sin(q12(ii))+L3*sin(q123(ii));
                     0,            0,                                               1];
    H04 = [cos(q4(ii)), -sin(q4(ii)), -Dx;
           sin(q4(ii)),  cos(q4(ii)),  Dy;
                     0,            0,   1];
    
    % Transform geometry
    DoorT  = H04*door;
    Link1T = H01*link1;
    Link2T = H02*link2;
    Link3T = H03*link3;
    
    % Plot geometry
    hold on
    if showpath == 1
        plot(x(1:ii),y(1:ii),'color',cmap(1,:),'linewidth',2)
        plot(x(1),y(1),'.','color',cmap(2,:),'markersize',18)
    end
    fill(Link3T(1,:),Link3T(2,:),col,'EdgeColor','k','linewidth',lw)
    fill(Link2T(1,:),Link2T(2,:),col,'EdgeColor','k','linewidth',lw)
    fill(Link1T(1,:),Link1T(2,:),col,'EdgeColor','k','linewidth',lw)
    fill(DoorT(1,:),DoorT(2,:),col,'EdgeColor','k','linewidth',lw)
    
    fill(hinge(1,:)-Dx,hinge(2,:)+Dy,col,'EdgeColor','k','linewidth',lw)
    fill(circ(1,:)/r*w4/2-Dx,circ(2,:)/r*w4/2+Dy,'k')
    fill(base(1,:),base(2,:),col,'EdgeColor','k','linewidth',lw)
    fill(circ(1,:)+H02(1,3),circ(2,:)+H02(2,3),'k')
    fill(circ(1,:)+H01(1,3),circ(2,:)+H01(2,3),'k')
    fill(circ(1,:),circ(2,:),'k')
    
    plot([-7 7]*dx,-2*1.2*r*[1 1],'k','linewidth',lw)
    for jj = -7:7-1
        plot(([0 1]+jj)*dx,[-2*1.2*r-dx,-2*1.2*r],'k','linewidth',lw)
    end
    
    wall1 = fill([-0.2 min(hinge(1,:)) min(hinge(1,:)) -0.2 -0.2]-Dx,3/4*[w4 w4 -w4 -w4 w4]+Dy,'w','EdgeColor','k','linewidth',lw);
    wall2 = fill(L4+[0.2 0.01 0.01 0.2 0.2]-Dx,3/4*[w4 w4 -w4 -w4 w4]+Dy,'w','EdgeColor','k','linewidth',lw);
    
    if showpath == 1
        plot(x(ii),y(ii),'.','color',cmap(1,:),'markersize',22)
    end
    
    % Set axis and figure properties
    axis('equal');
    axis([-0.1-Dx, L4-Dx+0.1, -0.2, 0.8])
    xlabel('x')
    ylabel('y')
    
    set(gca,'Nextplot','Replacechildren')
    set(gcf,'Position',[1300 550 592 564])
    hatch(wall1,[45,5,lw],'k')
    hatch(wall2,[45,5,lw],'k')
    
    % Capture frame
    if rec == 1
        writeVideo(writerObj,getframe(gcf));
    else
        F{ii} = getframe;
    end
end

if rec == 1
    close(writerObj)
end