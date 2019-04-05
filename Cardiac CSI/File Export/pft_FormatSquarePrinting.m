function pft_FormatSquarePrinting(HandleToFigure, HandleToAxes)

% Set the size of the figure and paper for output, regardless of the screen size
ScreenSize = get(0, 'ScreenSize');

wd = ScreenSize(3);
ht = ScreenSize(4);

Position = [ 0.1*wd, 0.1*ht, 0.8*wd, 0.8*ht ];

A4HT = 29.7;        % In cm
A4WD = 21.0;        % In cm
Edge = 1.00;        % In cm
Diff = A4HT - A4WD;

set(HandleToFigure, 'NextPlot', 'add');
set(HandleToFigure, 'PaperPositionMode', 'manual');
set(HandleToFigure, 'PaperType', 'a4');
set(HandleToFigure, 'PaperOrientation', 'landscape');
set(HandleToFigure, 'Position', Position);
set(HandleToFigure, 'PaperPosition', [Edge + Diff/2.0, Edge, A4WD - 2.0*Edge, A4WD - 2.0*Edge], 'PaperUnits', 'centimeters');

xt = get(HandleToAxes, 'XTick');
set(HandleToAxes, 'FontSize', 16);

yt = get(HandleToAxes, 'YTick');
set(HandleToAxes, 'FontSize', 16);

set(HandleToAxes, 'LineWidth', 1);

set(HandleToFigure, 'Color', [1 1 1]);
set(HandleToAxes, 'Color', [1 1 1]);

end

