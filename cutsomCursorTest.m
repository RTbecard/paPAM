t=0:.01:7; hp=plot(t,sin(t));
hCursorbar = graphics.cursorbar(hp); drawnow
hCursorbar.CursorLineColor = [.9,.3,.6]; % default=[0,0,0]='k'
hCursorbar.CursorLineStyle = '-';       % default='-'
hCursorbar.CursorLineWidth = 2.5;        % default=1
hCursorbar.Orientation = 'vertical';     % =default
hCursorbar.TargetMarkerSize = 12;        % default=8
hCursorbar.TargetMarkerStyle = 'x';      % default='s' (square)
