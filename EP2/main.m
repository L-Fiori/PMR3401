close all; clear all; clc

x = 22;
y = 10;

dx = 1;
dy = 1;

Az = zeros((y/dy) + 1, (x/dx) + 1);

Az_new = itemA(Az, dx, dy);

h = heatmap(Az_new, 'xlabel', 'x', 'ylabel', 'y', 'Colormap', hot);

h.GridVisible = 'off';
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
