function tue = tue_color(demo)
% TUE_COLOR defines a struct with colors, used by the Eindhoven University
% of Technology (TU/e).
% 
% TUE_COLOR('demo') also provides a demo plot with the possible colors.
%
% Possible colors:
%       r : red
%       db: dark blue
%       b : blue
%       c : cyan
%       g : green
%       y : yellow
%
% Example for using TU/e red in your own plot:
%       tue = TUE_COLOR;
%       plot(x,y,'Color',tue.r);
%
% (C) 2018 Stefan Driessen
 
switch nargin
    case 0
        d = false;
    case 1
        if strcmp(demo,'demo')
            d = true;
        else
            sprintf("Invalid argument %s",demo);
        end
end

tue = struct('r' , [0.784 0.098 0.098],...
             'db', [0.063 0.063 0.451],...
             'b' , [0     0.4   0.8  ],...
             'c' , [0     0.635 0.871],...
             'g' , [0.518 0.824 0    ],...
             'y' , [0.808 0.875 0    ]);

if d==true
    x = 0:0.001:4*pi;
    y = sin(x);
    step = 0.5;
    tue_cell = struct2cell(tue);
    figure('Name','TU/e color plot demo');
    hold on;
    for i=1:size(tue_cell,1)
        plot(x+(i-1)*step,y,'Color',tue_cell{i});
    end
    legend('red (r)','dark blue (db)','blue (b)','cyan (c)','green (g)','yellow (y)')
end
