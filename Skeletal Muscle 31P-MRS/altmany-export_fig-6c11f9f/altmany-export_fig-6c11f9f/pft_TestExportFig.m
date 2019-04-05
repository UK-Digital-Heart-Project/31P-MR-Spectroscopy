clear all
close all
clc

fclose('all');

Image = imread('pout.tif');

Dims = size(Image);
Rows = Dims(1);
Cols = Dims(2);

Width  = Cols;
Height = Rows;

iptsetpref('ImshowBorder', 'tight');

f = figure('Name', 'Testing Axes Properties', 'MenuBar', 'none', 'NumberTitle', 'off', 'Units', 'pixels');

imshow(Image, [ min(Image(:)), max(Image(:))]);

a = gca;

% colorbar('peer', a, 'Location', 'EastOutside');

% In the axes position [1, 1, 240, 291], the origin is 1-based and the width and height are inclusive pixel counts
set(a, 'Units', 'pixels');
AxesPosition = get(a, 'Position');

% The image is stored at its correct native size with a tight axes border within the main figure
% When the border is loose, there is an extra pixel in each dimension
ExportFileName = 'C:\Users\ptokarcz\Desktop\TestExportFig\Test-No-Colorbar';
export_fig(a, ExportFileName, '-png', '-jpg', '-bmp', '-tif');

% The image is stored at its correct native size with a tight axes border within the main figure
% When the border is loose, there is an extra pixel in each dimension
g = getframe(a);
[ Image, Map ] = frame2im(g);
imwrite(Image, 'C:\Users\ptokarcz\Desktop\TestExportFig\Test-No-Colorbar-Manual-Test.png');
imwrite(Image, 'C:\Users\ptokarcz\Desktop\TestExportFig\Test-No-Colorbar-Manual-Test.jpg');
imwrite(Image, 'C:\Users\ptokarcz\Desktop\TestExportFig\Test-No-Colorbar-Manual-Test.bmp');
imwrite(Image, 'C:\Users\ptokarcz\Desktop\TestExportFig\Test-No-Colorbar-Manual-Test.tif');

click = waitforbuttonpress;

iptsetpref('ImshowBorder', 'loose');