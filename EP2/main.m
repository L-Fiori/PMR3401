close all; clear all; clc

% Itens a) e b)

% Tamanho da malha utilizada (metade do dominio original,
% pois o problema é simétrico)
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
plotA(Az, title_def, M, N, dx)

% Item c)

[X, Y, Z] = meshgrid(1:N, 1:2*M-1, 1);

% Montagem da matriz de mis

mi = zeros(2*M-1, N);

mi_0 = 4*pi*10^-7;
mi_ferro = 2500*mi_0;
mi_ar = 1*mi_0;
mi_b = mi_ar;

for j=2:2*M-2
    for i=2:(N-1)
        if ((i-1)*dx =< 4)
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

[B_x, B_y, H_x, H_y] = itemC(Az, B_x, B_y, dx, dy, N, mi);

step = 1/dx;

titledef_B = 'Campo vetorial B (T) - passo ';
titledef_H = 'Campo vetorial H (A/m^2) - passo ';

plotC(X, Y, B_x, B_y, H_x, H_y, step, titledef_B, titledef_H)

% Item d)

f_x = zeros(length(Az));
f_y = zeros(length(Az));

n = length(Az);

[F_x, F_y] = itemD(B_x, B_y, dx, M, n, mi_0)  

% Item e1)

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
plotA(Az, title_def, M, N, dx)

% "item c" do e1
B_x = zeros(2*M - 1, N);
B_y = zeros(2*M - 1, N);

[B_x, B_y, H_x, H_y] = itemC(Az, B_x, B_y, dx, dy, N, mi);

step = 1/dx;

plotC(X, Y, B_x, B_y, H_x, H_y, step)

% "item d" do e1

f_x = zeros(length(Az));
f_y = zeros(length(Az));

[F_x, F_y] = itemD(B_x, B_y, dx, M, n, mi_0)

% item e2)

% "itens a e b" do e2


