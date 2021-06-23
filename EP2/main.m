close all; clear all; clc

x = 22;
y = 10;

dx = 1;
dy = 1;

M = (y/dy) + 1;
N = (x/dx) + 1;

Az = zeros(M, N);

Az_new = itemA(Az, dx, dy);

Az_full = zeros(2*M, N);

Az_full(M+1:2*M, :) = Az_new(:, :);

for j=1:M
    for i=(2:N-1)
        Az_full(j, i) = Az_new(M+1-j, i);
    end
end 


h = heatmap(Az_full, 'xlabel', 'x', 'ylabel', 'y', 'Colormap', flipud(hot));

h.GridVisible = 'off';
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));


% Item c)



