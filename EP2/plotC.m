function plotC(X, Y, B_x, B_y, H_x, H_y, step, titledef_B, titledef_H)

figure
q = quiver(X(1:step:end, 1:step:end), Y(1:step:end, 1:step:end), B_x(1:step:end, 1:step:end), B_y(1:step:end, 1:step:end));
q.AutoScaleFactor = 1.4;

title([titledef_B, num2str(dx)],'FontSize', 12);
set(gca,'XTickLabel',{0, 4, 8, 12, 16, 20})
set(gca,'YTickLabel',{-8, -6, -4, -2, 0, 2, 4, 6, 8, 10})
xticks(linspace(0, N, 12));
xticklabels(arrayfun(@num2str, (linspace(0, 22, 12)), 'UniformOutput', false));
yticks(linspace(0, 2*M-1, 11));
yticklabels(arrayfun(@num2str, (linspace(-10, 10, 11)), 'UniformOutput', false));
ylabel('y (cm)') 
xlabel('x (cm)')
axis equal

figure
q_h = quiver(X(1:step:end, 1:step:end), Y(1:step:end, 1:step:end), H_x(1:step:end, 1:step:end), H_y(1:step:end, 1:step:end));
q_h.AutoScaleFactor = 3;

title([titledef_H, num2str(dx)],'FontSize', 12);
set(gca,'XTickLabel',{0, 4, 8, 12, 16, 20})
set(gca,'YTickLabel',{-8, -6, -4, -2, 0, 2, 4, 6, 8, 10})
xticks(linspace(0, N, 12));
xticklabels(arrayfun(@num2str, (linspace(0, 22, 12)), 'UniformOutput', false));
yticks(linspace(0, 2*M-1, 11));
yticklabels(arrayfun(@num2str, (linspace(-10, 10, 11)), 'UniformOutput', false));
ylabel('y (cm)')
xlabel('x (cm)')
axis equal


end