close all; clear all; clc

% Itens a) e b)

% Tamanho da malha utilizada (metade do dominio original,
% pois o problema nao eh simetrico)
x = 22;
y = 10;

% Especificacao do passo
dx = 0.1;
dy = 0.1;

M = (y/dy) + 1;
N = (x/dx) + 1;

Az = zeros(M, N);

Az_new = itemA(Az, dx, dy);

Az = zeros(2*M-1, N);

Az(M:end, :) = Az_new(:, :);

for j=1:M-1
    for i=(2:N-1)
        Az(j, i) = Az_new(M+1-j, i);
    end
end 

title_def = 'a) Az (Wb/m) - passo ';
isE2 = false;
isE = false;
plotA(Az, title_def, M, N, dx, isE2)

% Item c)

[X, Y, Z] = meshgrid(1:N, 1:2*M-1, 1);
x = 1:N;
y = 1:(2*M-1);


% Montagem da matriz de mis

mi = zeros(2*M-1, N);

mi_0 = 4*pi*10^-7;
mi_ferro = 2500*mi_0;
mi_ar = 1*mi_0;
mi_b = mi_ar;

mi_ferro_x = 1200*mi_0;
mi_ferro_y = 2500*mi_0;

for j=2:2*M-2
    for i=2:(N-1)
        if ((i-1)*dx <= 4)
            mi(j, i) = 1/mi_ferro;
        elseif((i-1)*dx > 4 && (i-1)*dx < 5)
            mi(j, i) = 1/mi_ar;
        elseif((i-1)*dx > 5 && (i-1)*dx < 14)
            if ((j-1)*dy < 4 || (j-1)*dy > 16)
                mi(j, i) = 1/mi_ferro;
            elseif((j-1)*dy > 4 && (j-1)*dy < 16)
                mi(j, i) = 1/mi_ar;
            end
        elseif((i-1)*dx > 14 && (i-1)*dx < 16)
            if ((j-1)*dy < 4 || (j-1)*dy > 16)
                mi(j, i) = 1/mi_ferro;
            elseif((j-1)*dy > 4 && (j-1)*dy < 16)
                mi(j, i) = 1/mi_b;
            end
        elseif((i-1)*dx > 16 && (i-1)*dx < 20)
            mi(j, i) = 1/mi_ferro;
        elseif((i-1)*dx > 20 && (i-1)*dx < 22)
            if((j-1)*dy > 4 && (j-1)*dy < 16)
                mi(j, i) = 1/mi_b;
            end
        end
    end
end

B_x = zeros(2*M - 1, N);
B_y = zeros(2*M - 1, N);

[B_x, B_y, H_x, H_y] = itemC(Az, B_x, B_y, dx, dy, N, mi, mi, isE)

step = 1/dx;

titledef_B = 'Campo vetorial B (T) - passo ';
titledef_H = 'Campo vetorial H (A/m^2) - passo ';

plotC(X, Y, B_x, B_y, H_x, H_y, step, titledef_B, titledef_H, dx, M, N)

% Item d)

f_x = zeros(2*M-1);
f_y = zeros(2*M-1);

n = 2*M-1;

[F_x, F_y] = itemD(B_x, B_y, dx, M, n, mi_0)

% Item e1)
isE2 = false;
isE = true;

% "itens a e b" do e1
Az = zeros(M, N);
Az_new = itemE1(Az, dx, dy);

Az = zeros(2*M-1, N);

Az(M:end, :) = Az_new(:, :);

for j=1:M-1
    for i=(2:N-1)
        Az(j, i) = Az_new(M+1-j, i);
    end
end

title_def = 'e1) Az (Wb/m) - passo ';
plotA(Az, title_def, M, N, dx, isE2)

% "item c" do e1
B_x = zeros(2*M - 1, N);
B_y = zeros(2*M - 1, N);

mi_x = zeros(2*M-1, N);
mi_y = mi;

for j=2:2*M-2
    for i=2:(N-1)
        if ((i-1)*dx <= 4)
            mi_x(j, i) = 1/mi_ferro_x;
        elseif((i-1)*dx > 4 && (i-1)*dx < 5)
            mi_x(j, i) = 1/mi_ar;
        elseif((i-1)*dx > 5 && (i-1)*dx < 14)
            if ((j-1)*dy < 4 || (j-1)*dy > 16)
                mi_x(j, i) = 1/mi_ferro_x;
            elseif((j-1)*dy > 4 && (j-1)*dy < 16)
                mi_x(j, i) = 1/mi_ar;
            end
        elseif((i-1)*dx > 14 && (i-1)*dx < 16)
            if ((j-1)*dy < 4 || (j-1)*dy > 16)
                mi_x(j, i) = 1/mi_ferro_x;
            elseif((j-1)*dy > 4 && (j-1)*dy < 16)
                mi_x(j, i) = 1/mi_b;
            end
        elseif((i-1)*dx > 16 && (i-1)*dx < 20)
            mi_x(j, i) = 1/mi_ferro_x;
        elseif((i-1)*dx > 20 && (i-1)*dx < 22)
            if((j-1)*dy > 4 && (j-1)*dy < 16)
                mi_x(j, i) = 1/mi_b;
            end
        end
    end
end

[B_x, B_y, H_x, H_y] = itemC(Az, B_x, B_y, dx, dy, N, mi_x, mi_y, isE);

step = 1/dx;

plotC(X, Y, B_x, B_y, H_x, H_y, step, titledef_B, titledef_H, dx, M, N)

% "item d" do e1

f_x = zeros(length(Az));
f_y = zeros(length(Az));
n = 2*M - 1;

[F_x, F_y] = itemD(B_x, B_y, dx, M, n, mi_0)

% item e2)
isE2 = true;

% "itens a e b" do e2

dt = 0.000001; % passo = 0.1
t_range = 500; % Intervalo temporal para plotagem da força (*deltat) [500]
divs = 10; % Divisões do intervalo acima [10]
tamanho = t_range/divs;

Az = zeros(M, N, 2); % 50 ~~~~~~
Az_new = itemE2(Az, dx, dt, t_range, divs, M, N);
    
Az = zeros(2*M-1, N, 3);

Az(M:end, :, 1) = Az_new(:, :, 1);
Az(M:end, :, 2) = Az_new(:, :, 10);
Az(M:end, :, 3) = Az_new(:, :, 50);

for j=1:M-1
    for i=(2:N-1)
        Az(j, i, 1) = Az_new(M+1-j, i, 1);
        Az(j, i, 2) = Az_new(M+1-j, i, 10);
        Az(j, i, 3) = Az_new(M+1-j, i, 50);
    end
end 

titledef = 'e2) Az (Wb/m) - passo ';
isE2 = true;
num = 10;
plotA(Az(:, :, 1), titledef, M, N, dx, isE2, dt, num)

num = 100;
plotA(Az(:, :, 2), titledef, M, N, dx, isE2, dt, num)

num = 500;
plotA(Az(:, :, 3), titledef, M, N, dx, isE2, dt, num)

% "item c" do item e2
B_x = zeros(2*M - 1, N);
B_y = zeros(2*M - 1, N);

isE = true;
[B_x, B_y, H_x, H_y] = itemC(Az(:, :, 1), B_x, B_y, dx, dy, N, mi_x, mi_y, isE);

titledef_B = 'Campo vetorial B (T) - passo ';
titledef_H = 'Campo vetorial H (A/m^2) - passo ';
step = 1/dx-4;
plotC(X, Y, B_x, B_y, H_x, H_y, step, titledef_B, titledef_H, dx, M, N)

B_x = zeros(2*M - 1, N);
B_y = zeros(2*M - 1, N);

[B_x, B_y, H_x, H_y] = itemC(Az(:, :, 2), B_x, B_y, dx, dy, N, mi_x, mi_y, isE);

plotC(X, Y, B_x, B_y, H_x, H_y, step, titledef_B, titledef_H, dx, M, N)

B_x = zeros(2*M - 1, N);
B_y = zeros(2*M - 1, N);

[B_x, B_y, H_x, H_y] = itemC(Az(:, :, 3), B_x, B_y, dx, dy, N, mi_x, mi_y, isE);

plotC(X, Y, B_x, B_y, H_x, H_y, step, titledef_B, titledef_H, dx, M, N)

% "item d" do item e2

Az = zeros(2*M-1, N, 50);

Az(M:end, :, :) = Az_new(:, :, :);

for j=1:M-1
    for i=(2:N-1)
        Az(j, i, :) = Az_new(M+1-j, i, :);
    end
end 

B_x = zeros(2*M - 1, N, 50);
B_y = zeros(2*M - 1, N, 50);


for k=1:tamanho
    [B_x(k), B_y(k), H_x(k), H_y(k)] = itemC(Az(:, :, k), B_x(:, :, k), B_y(:, :, k), dx, dy, N, mi_x, mi_y, isE);
end

n = 2*M-1;
F_x = zeros(n, 50);
F_y = zeros(n, 50);

for k=1:tamanho
    [F_x(k), F_y(k)] = itemD(B_x(k), B_y(k), dx, M, n, mi_0); 
end

t = linspace(divs*dt, t_range*dt, tamanho);

figure
plot(t, F_x);
title(['Forca horizontal aplicada na armadura, passo=', num2str(dx),'cm, \Deltat=', num2str(dt),'s']);
xlabel('Tempo (s)');
ylabel('Fx (N)');
grid on;