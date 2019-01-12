function PlotHybTrajSnap(traj,sgnls,grd,tcut,options)

% The function PlotHybTrajSnap plots the supplied hybrid trajectories.
%
% PlotHybTraj(traj,sgnls,grd) plots for the hybrid trajectory/trajectories 
% traj (which should be a cell array if multiple trajectories should be 
% plotted) the signals sgnls.x and sgnls.u (indicating the states and
% inputs to be shown). The signals are depicted in the subplots in a grid
% grd, e.g. grd = [2,2] will result in a 2 x 2 subplot figure. The inputs
% are plotted first, followed by the states. If the signals should be shown
% in a different order, one can supply a structure grd with grd.x and grd.u
% matrices where each row contains the location of the state (resp. input)
% that should be plotted in the grid of subplots, e.g. grd.x = [1,2]
% indicates that the state is plotted in the (1,2) element of subplots (and
% that only one state is plotted, the state in sgnls.x).
%
% PlotHybTraj(traj,sgnls,grd,options) plots for the hybrid trajectories
% with selected plotting options. The options can be specified by supplying 
% a struct (options) with fields corresponding to the options that are 
% desired. The options that can be specified are:
%
% barcolor        :  A matrix specifying the colors corresponding to
%                    different modes in constructing the event counter bar. 
%                    The (i)th row contains the RGB color parameters 
%                    corresponding to mode i. 
% barspos         :  String indicating if the event counter bar(s) should 
%                    be plotted above ('top') or below ('bottom') the 
%                    trajectory plots.
% cntr            :  Illustrate the counter on the counter bar(s) by means
%                    of numbers ('num') or changing color shades ('color')
% cntrbar         :  Vector with 1's and 0's that specify if for the
%                    corresponding trajectory an event counter bar should 
%                    be plotted yes or no, respectively.
% dotsu           :  Add dots to the input signals at the event times 
%                    yes (1) or no (0). 
% evntlines       :  Add vertical lines to the plots at the event times 
%                    yes (1) or no (0).
% font            :  Font used for the legend, axis labels, counters, etc.
% fontsize.axes   :  Font size (in pt's) of the axes annotations (i.e. the 
%                    X and Y Tick labels).
% fontsize.labels :  Font size (in pt's) of the axes labels.
% fontsize.text   :  Font size (in pt's) of the event counters and legend. 
% grid            :  Add a grid to the different trajectory state and input 
%                    plots minor (2) normal (1) or no (0).
% Hbar            :  Height (relative, i.e. between 0 and 1) of the
%                    horizontal bars illustrating the event counter as a 
%                    function of time.
% height          :  Heigth of the figure in centimeters.
% labels          :  Cell array (1 x N) containing the ylabels of the 
%                    different subplots ordered along the rows of the 
%                    subplot (e.g. the 3rd label will be added to the (2,1) 
%                    entry in a 2x2 subplot structure.
% legend          :  Cell array containing the legend entries corresponding 
%                    to the different trajectories.
% linecolor       :  Cell array with linecolor identifiers (either 1 x 3 
%                    RGB vector or string as 'b' for example) corresponding
%                    to the different trajectories.
% linestyle       :  Cell array with linestyle identifiers (e.g. '-.')
%                    corresponding to the different trajectories.
% linewidth       :  Vector with the linewidth of the different 
%                    trajectories plotted.    
% marginx         :  Margin (relative, i.e. between 0 and 1) around the 
%                    figure (subplots) in x-direction.    
% marginy         :  Margin (relative, i.e. between 0 and 1) around the 
%                    figure (subplots) in y-direction.    
% marker          :  String identifier for the type of marker used for 
%                    illustrating the left and right limits at the event 
%                    times of the hybrid trajectories. See the help of 
%                    plot.m for possible markers.
% markersize      :  The size of the markers illustrating the state and 
%                    input at the event times.
% showlegend      :  Add a legend to the figure yes (1) or no (0). The 
%                    legend will be added to the first trajectory subplot.
% spacingx        :  Spacing (relative, i.e. between 0 and 1) between two
%                    consecutive subplots in x-direction.
% spacingy        :  Spacing (relative, i.e. between 0 and 1) between two
%                    consecutive subplots in y-direction. This is also the
%                    scaled distance between the counter bar plots and the
%                    trajectory subplots.
% spacingb        :  Spacing (relative, i.e. between 0 and 1) between two
%                    consecutive event counter bars in y-direction.
% width           :  Width of the figure in centimeters.
%
% See also HybridSim, PlotHybTrajOptns (the default settings)
% M.W.L.M. Rijnen 
% December 12, 2016 
% Eindhoven University of Technology

% Determine ordering of the subplots
if isfield(grd,'x') 
    if isfield(grd,'u')
        N_rows = max([grd.u;grd.x]*[1;0]);
        N_cols = max([grd.u;grd.x]*[0;1]);
        ID_x   = zeros(N_rows,N_cols);
        ID_u   = zeros(N_rows,N_cols);
        for i = 1:size(grd.x,1)
            ID_x(grd.x(i,1),grd.x(i,2)) = sgnls.x(i);
        end
        for i = 1:size(grd.u,1)
            ID_u(grd.u(i,1),grd.u(i,2)) = sgnls.u(i);
        end
    else
        error('If the plot grid for the states is supplied it should also be supplied for the inputs u.')
    end
else
    N_rows = grd(1);
    N_cols = grd(2);
    ID_x   = zeros(N_cols,N_rows);
    ID_u   = zeros(N_cols,N_rows);
    for i = 1:length(sgnls.u)
        ID_u(i) = sgnls.u(i);
    end
    for i = 1:length(sgnls.x)
        ID_x(length(sgnls.u)+i) = sgnls.x(i);
    end
    ID_x = ID_x.';
    ID_u = ID_u.';
end

% Set the plot options
if ~iscell(traj)
    traj = {traj};
end
Nm = 1;
for i = 1:length(traj)
    if max(traj{i}.m) > Nm
        Nm = max(traj{i}.m);
    end
end
if nargin == 4
    options = PlotHybTrajOptns(length(traj),Nm,ID_x,ID_u);
elseif nargin == 5
    options = PlotHybTrajOptns(length(traj),Nm,ID_x,ID_u,options);
else
    error('Incorrect number of input arguments to PlotHybTraj.m.')
end
i_bars = find(options.cntrbar == 1);
N_bars = length(i_bars);

% Apply the font settings
set(0,'DefaultAxesFontName', options.font)
set(0,'DefaultTextFontName', options.font)
set(0,'DefaultAxesFontSize', options.fontsize.axes)
set(0,'DefaultTextFontSize', options.fontsize.text)

% Dimensions and position of the subplots
spacingx = options.spacingx;
spacingy = options.spacingy;
spacingb = options.spacingb;
marginx  = options.marginx;
marginy  = options.marginy;
Hbar     = options.Hbar;
W        = (1-1.2*marginx-(N_cols-1)*spacingx)/N_cols;
H        = (1-1.5*marginy-N_bars*Hbar-(N_bars>0)*(N_bars-1)*spacingb-(N_rows-1+(N_bars>0))*spacingy)/N_rows;


% Plot hybrid trajectory/trajectories ************************************%

% Plot counter bar(s)
h = zeros((N_rows+N_bars)*N_cols,1);
ind = 0;
hold on
for i = 1:N_bars
    if isfield(traj{i_bars(i)},'te')
        t_vec = [traj{i_bars(i)}.t{1}(1);traj{i_bars(i)}.te;traj{i_bars(i)}.t{end}(end)];
    elseif length(traj{i_bars(i)}.t) > 1
        error('Event times are not included in the hybrid trajectory event though there are multiple trajectory segments.')
    else
        t_vec = [traj{i_bars(i)}.t{1}(1);traj{i_bars(i)}.t{end}(end)];
    end
    for j = 1:N_cols
        ind    = (i-1)*N_cols+j; 
        h(ind) = subplot(N_rows+N_bars,N_cols,ind);
        plot([0 0],[-1 -1],'color',options.linecolor{i_bars(i)},...
                           'LineStyle',options.linestyle{i_bars(i)},...
                           'LineWidth',options.linewidth(i_bars(i)))
        hold on
        if strcmpi(options.cntr,'num')
            for k = 1:length(t_vec)-1
                fill([t_vec(k) t_vec(k+1) t_vec(k+1) t_vec(k)],[0 0 1 1],options.barcolor(traj{i_bars(i)}.m(k),:),'edgecolor','k')
                brightness = options.barcolor(traj{i_bars(i)}.m(k),:)*[299;587;114]/1000;
                if brightness < 0.5
                    text(0.5*(t_vec(k+1)+t_vec(k)),0.5,sprintf('%d',k-1),'VerticalAlignment','middle','HorizontalAlignment','center','color','w')
                else
                    text(0.5*(t_vec(k+1)+t_vec(k)),0.5,sprintf('%d',k-1),'VerticalAlignment','middle','HorizontalAlignment','center')
                end
                set(gca,'XTick',[])
            end
        elseif strcmpi(options.cntr,'color')        
            set(h(ind),'XColor','w','YColor','w');
            for k = 1:length(t_vec)-1
                fill([t_vec(k) t_vec(k+1) t_vec(k+1) t_vec(k)],[0 0 1 1],1+(options.barcolor(traj{i_bars(i)}.m(k),:)-1)*k/(length(t_vec)-1),'edgecolor','w','linewidth',1)
                set(gca,'XTick',[])
            end
        end
        ylim([0 1])
    end
end   
indb = ind;

% Plot trajectories
h_legend = zeros(1,length(traj));
for i = 1:N_rows
    for j = 1:N_cols
        ind    = indb + (i-1)*N_cols+j; 
        h(ind) = subplot(N_rows+N_bars,N_cols,ind);
        hold on
        
        if ID_x(i,j) ~= 0 && ID_u(i,j) == 0
            ID = ID_x(i,j);
        
            for k = 1:length(traj)
                for ii = 1:length(traj{k}.t)
                    if ID <= size(traj{k}.x{ii},2)-1
                        if traj{k}.t{ii}(end) <= tcut
                            if i == 1 && j == 1 && ii == 1
                                h_legend(k) = plot(traj{k}.t{ii},traj{k}.x{ii}(:,ID),...
                                                   'color',options.linecolor{k},...
                                                   'linestyle',options.linestyle{k},...
                                                   'linewidth',options.linewidth(k));
                            else
                                plot(traj{k}.t{ii},traj{k}.x{ii}(:,ID),...
                                     'color',options.linecolor{k},...
                                     'linestyle',options.linestyle{k},...
                                     'linewidth',options.linewidth(k))
                            end
                            if ii == 1 && length(traj{k}.t) > 1
                                plot(traj{k}.te(ii),traj{k}.xel{ii}(ID),...
                                     'color',options.linecolor{k},...
                                     'marker',options.marker{k,1},...
                                     'linestyle','none',...
                                     'markersize',options.markersize(k,1))
                            elseif ii < length(traj{k}.t)
                                plot(traj{k}.te(ii),traj{k}.xel{ii}(ID),...
                                     'color',options.linecolor{k},...
                                     'marker',options.marker{k,1},...
                                     'linestyle','none',...
                                     'markersize',options.markersize(k,1))
                                plot(traj{k}.te(ii-1),traj{k}.xer{ii-1}(ID),...
                                     'color',options.linecolor{k},...
                                     'marker',options.marker{k,2},...
                                     'linestyle','none',...
                                     'markersize',options.markersize(k,2))
                            elseif length(traj{k}.t) > 1
                                plot(traj{k}.te(ii-1),traj{k}.xer{ii-1}(ID),...
                                     'color',options.linecolor{k},...
                                     'marker',options.marker{k,2},...
                                     'linestyle','none',...
                                     'markersize',options.markersize(k,2))
                            end
                        else
                            jj = find(traj{k}.t{ii}>tcut,1,'first');
                            if i == 1 && j == 1 && ii == 1
                                h_legend(k) = plot(traj{k}.t{ii}(1:jj),traj{k}.x{ii}(1:jj,ID),...
                                                   'color',options.linecolor{k},...
                                                   'linestyle',options.linestyle{k},...
                                                   'linewidth',options.linewidth(k));
                            else
                                plot(traj{k}.t{ii}(1:jj),traj{k}.x{ii}(1:jj,ID),...
                                     'color',options.linecolor{k},...
                                     'linestyle',options.linestyle{k},...
                                     'linewidth',options.linewidth(k))
                            end
                            if ii == 1 && length(traj{k}.t) > 1
                                plot(traj{k}.te(ii),traj{k}.xel{ii}(ID),...
                                     'color',options.linecolor{k},...
                                     'marker',options.marker{k,1},...
                                     'linestyle','none',...
                                     'markersize',options.markersize(k,1))
                            elseif ii < length(traj{k}.t)
                                plot(traj{k}.te(ii),traj{k}.xel{ii}(ID),...
                                     'color',options.linecolor{k},...
                                     'marker',options.marker{k,1},...
                                     'linestyle','none',...
                                     'markersize',options.markersize(k,1))
                                plot(traj{k}.te(ii-1),traj{k}.xer{ii-1}(ID),...
                                     'color',options.linecolor{k},...
                                     'marker',options.marker{k,2},...
                                     'linestyle','none',...
                                     'markersize',options.markersize(k,2))
                            elseif length(traj{k}.t) > 1
                                plot(traj{k}.te(ii-1),traj{k}.xer{ii-1}(ID),...
                                     'color',options.linecolor{k},...
                                     'marker',options.marker{k,2},...
                                     'linestyle','none',...
                                     'markersize',options.markersize(k,2))
                            end
                        end
                    end
                end           
            end
            
        elseif ID_u(i,j) ~= 0 && ID_x(i,j) == 0
            
            ID = ID_u(i,j);
            for k = 1:length(traj)
                for ii = 1:length(traj{k}.t)
                    if ID <= size(traj{k}.u{ii},2)
                        if traj{k}.t{ii}(end) <= tcut
                            if i == 1 && j == 1 && ii == 1
                                h_legend(k) = plot(traj{k}.t{ii},traj{k}.u{ii}(:,ID),...
                                                   'color',options.linecolor{k},...
                                                   'linestyle',options.linestyle{k},...
                                                   'linewidth',options.linewidth(k));
                            else
                                plot(traj{k}.t{ii},traj{k}.u{ii}(:,ID),...
                                     'color',options.linecolor{k},...
                                     'linestyle',options.linestyle{k},...
                                     'linewidth',options.linewidth(k))
                            end
                            if options.dotsu == 1
                                if ii == 1 && length(traj{k}.t) > 1
                                    plot(traj{k}.te(ii),traj{k}.u{ii}(traj{k}.t{ii}==traj{k}.te(ii),ID),...
                                         'color',options.linecolor{k},...
                                         'marker',options.marker{k,1},...
                                         'linestyle','none',...
                                         'markersize',options.markersize(k,1))
                                elseif ii < length(traj{k}.t)
                                    plot(traj{k}.te(ii),traj{k}.u{ii}(traj{k}.t{ii}==traj{k}.te(ii),ID),...
                                         'color',options.linecolor{k},...
                                         'marker',options.marker{k,1},...
                                         'linestyle','none',...
                                         'markersize',options.markersize(k,1))
                                    plot(traj{k}.te(ii-1),traj{k}.u{ii}(traj{k}.t{ii}==traj{k}.te(ii-1),ID),...
                                         'color',options.linecolor{k},...
                                         'marker',options.marker{k,2},...
                                         'linestyle','none',...
                                         'markersize',options.markersize(k,2))
                                elseif length(traj{k}.t) > 1
                                    plot(traj{k}.te(ii-1),traj{k}.u{ii}(traj{k}.t{ii}==traj{k}.te(ii-1),ID),...
                                         'color',options.linecolor{k},...
                                         'marker',options.marker{k,2},...
                                         'linestyle','none',...
                                         'markersize',options.markersize(k,2))
                                end
                            end
                        else
                            jj = find(traj{k}.t{ii}>tcut,1,'first');
                            if i == 1 && j == 1 && ii == 1
                                h_legend(k) = plot(traj{k}.t{ii}(1:jj),traj{k}.u{ii}(1:jj,ID),...
                                                   'color',options.linecolor{k},...
                                                   'linestyle',options.linestyle{k},...
                                                   'linewidth',options.linewidth(k));
                            else
                                plot(traj{k}.t{ii}(1:jj),traj{k}.u{ii}(1:jj,ID),...
                                     'color',options.linecolor{k},...
                                     'linestyle',options.linestyle{k},...
                                     'linewidth',options.linewidth(k))
                            end
                            if options.dotsu == 1
                                if ii == 1 && length(traj{k}.t) > 1
                                    plot(traj{k}.te(ii),traj{k}.u{ii}(traj{k}.t{ii}==traj{k}.te(ii),ID),...
                                         'color',options.linecolor{k},...
                                         'marker',options.marker{k,1},...
                                         'linestyle','none',...
                                         'markersize',options.markersize(k,1))
                                elseif ii < length(traj{k}.t)
                                    plot(traj{k}.te(ii),traj{k}.u{ii}(traj{k}.t{ii}==traj{k}.te(ii),ID),...
                                         'color',options.linecolor{k},...
                                         'marker',options.marker{k,1},...
                                         'linestyle','none',...
                                         'markersize',options.markersize(k,1))
                                    plot(traj{k}.te(ii-1),traj{k}.u{ii}(traj{k}.t{ii}==traj{k}.te(ii-1),ID),...
                                         'color',options.linecolor{k},...
                                         'marker',options.marker{k,2},...
                                         'linestyle','none',...
                                         'markersize',options.markersize(k,2))
                                elseif length(traj{k}.t) > 1
                                    plot(traj{k}.te(ii-1),traj{k}.u{ii}(traj{k}.t{ii}==traj{k}.te(ii-1),ID),...
                                         'color',options.linecolor{k},...
                                         'marker',options.marker{k,2},...
                                         'linestyle','none',...
                                         'markersize',options.markersize(k,2))
                                end
                            end
                        end
                    end
                end
            end
            
        elseif ID_u(i,j) ~= 0 && ID_x(i,j) ~= 0
            error('Both a state and an input are chosen to be plotted in the same subplot.')
        end
        
        if options.evntlines == 1
            ylims = get(h(ind),'YLim');
            for k = 1:length(traj)
                if isfield(traj{k},'te')
                    for ii = 1:length(traj{k}.te)
                        plot([traj{k}.te(ii) traj{k}.te(ii)],ylims,...
                                     'color',options.linecolor{k},...
                                     'linestyle',':',...
                                     'linewidth',0.5)
                    end
                end
            end
        end
               
        if options.grid == 1
            grid on
        elseif options.grid == 2
            grid minor 
        end
        
        ylabel(options.labels{ind-indb},'FontSize',options.fontsize.labels,'Interpreter','latex')
    end
end 
        
% Change appearance of the subplots
tmin = 1e10;
tmax = -1e10;
for i = 1:length(traj)
    if tmin > traj{i}.t{1}(1)
        tmin = traj{i}.t{1}(1);
    end
    if tmax < traj{i}.t{end}(end)
        tmax = traj{i}.t{end}(end);
    end
end
linkaxes(h,'x')
xlim([tmin tmax]);

for i = 1:N_bars
    for j = 1:N_cols
        ind    = (i-1)*N_cols+j; 
        if strcmpi(options.barspos,'top')
            if j == 1 && length(options.cntrbar) > 1
                legend(h(ind),'','Location',[2*marginx/3 1-marginy/2-i*Hbar-(i-1)*spacingb marginx/6 Hbar])
                legend(h(ind),'boxoff')
            end 
            set(h(ind),'Position',[marginx+(j-1)*(W+spacingx) 1-marginy/2-i*Hbar-(i-1)*spacingb W Hbar],...
                       'XTickLabel',[],'YTickLabel',[])
        elseif strcmpi(options.barspos,'bottom')
            if j == 1 && length(options.cntrbar) > 1
                legend(h(ind),'','Location',[2*marginx/3 1-marginy/2-N_rows*(H+spacingy)-i*Hbar-(i-1)*spacingb marginx/6 Hbar])
                legend(h(ind),'boxoff')
            end 
            set(h(ind),'Position',[marginx+(j-1)*(W+spacingx) 1-marginy/2-N_rows*(H+spacingy)-i*Hbar-(i-1)*spacingb W Hbar])
            if i < N_bars
                set(h(ind),'XTickLabel',[],'YTickLabel',[])
            else
                set(h(ind),'YTickLabel',[])
                axes(h(ind))
                xlabel('t [s]','FontSize',options.fontsize.labels,'Interpreter','latex')
            end
        end
        set(h(ind),'YTick',[])
        set(h(ind),'box','off','color','none')
        set(h(ind),'Layer','Top')
    end
end
for i = 1:N_rows
    for j = 1:N_cols
        ind    = indb + (i-1)*N_cols+j; 
        if strcmpi(options.barspos,'top')
            if i < N_rows
                set(h(ind),'XTicklabel',[])
            else
                xlabel('t [s]','FontSize',options.fontsize.labels,'Interpreter','latex')
            end
            set(h(ind),'Position',[marginx+(j-1)*(W+spacingx) 1-marginy/2-i*(H+spacingy)+(N_bars==0)*spacingy-N_bars*Hbar-(N_bars>0)*(N_bars-1)*spacingb W H])
        elseif strcmpi(options.barspos,'bottom')
            set(h(ind),'XTicklabel',[])
            set(h(ind),'Position',[marginx+(j-1)*(W+spacingx) 1-marginy/2-i*(H+spacingy)+spacingy W H])
        end
        if options.showlegend == 1 && i == 1 && j == 1
            legend(h_legend(1:length(options.legend)),options.legend,'FontSize',options.fontsize.text,'Interpreter','latex')
        end
    end
end

% Set figure dimensions
set(gca,'visible','off')
set(gcf,'units','centimeters')
set(gcf,'position',[4 2 options.width options.height])
end