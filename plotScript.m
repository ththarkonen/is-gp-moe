
clear
close all

[ fileName, folderPath] = uigetfile("*.mat");
filePath = fullfile( folderPath, fileName);

result = load( filePath );
result = result.ismoeObject;

xStar = result.X( 1, :);
yStar = result.Y( :, 1);

minY = min( yStar );
maxY = max( yStar );

[ nY, nX] = size( result.X );
[ highDensityMatrix, medianEstimate] = computeHighDensityArea( result, nX, nY);

figure();
hold on

hData = plot( result.data.x, result.data.y, '.');
hData.MarkerSize = 65;

h = imagesc( xStar, yStar, highDensityMatrix);

cmap = flipud( gray );
colormap( cmap );

h = plot( xStar, medianEstimate);
h.Color = [ 51, 51, 51] / 255;
h.LineWidth = 10;

hData = plot( result.data.x, result.data.y, '.', 'Color', hData.Color);
hData.MarkerSize = 65;

ylim([ minY, maxY]);

h = xlabel("$x$");
h.Interpreter = "latex";

h = ylabel("$y$");
h.Interpreter = "latex";

ax = gca();
ax.Layer = "top";
ax.YDir = "normal";
ax.FontSize = 55;
caxis([0, 4]);
