close all; clear all; clc

x = 22;
y = 10;

dx = 1;
dy = 1;

M = (y/dy) + 1;
N = (x/dx) + 1;

Az = zeros(M, N);

Az_new = itemA(Az, dx, dy);

Az = zeros(2*M, N);

Az(M+1:2*M, :) = Az_new(:, :);

for j=1:M
    for i=(2:N-1)
        Az(j, i) = Az_new(M+1-j, i);
    end
end 

%{
h = heatmap(Az, 'xlabel', 'x', 'ylabel', 'y', 'Colormap', flipud(hot));

h.GridVisible = 'off';
Ax = gca;
Ax.XDisplayLabels = nan(size(Ax.XDisplayData));
Ax.YDisplayLabels = nan(size(Ax.YDisplayData));
%}

% Item c)


%Implementacao utilizando curl()

[X, Y, Z] = meshgrid(1:N, 1:2*M-1, 1);

%{
X(:,:,2) = ones(2*M, N)+1;
Y(:,:,2) = ones(2*M, N)+1;
Z(:,:,2) = ones(2*M, N)+1;



F_x = zeros(2*M, N, 2);
F_y = zeros(2*M, N, 2);

Az(:, :, 2) = zeros(2*M, N);

[curl_x, curl_y, curl_z, cav] = curl(X, Y, Z, F_x, F_y, Az);

step = 1/dx;

q = quiver(X(1:step:end, 1:step:end), Y(1:step:end, 1:step:end), curl_x(1:step:end, 1:step:end), curl_y(1:step:end, 1:step:end));
q.AutoScaleFactor = 1.4;
%}

%Implementacao utilizando diferencas finitas

% Calculo de B

B_x = zeros(2*M - 1, N);
B_y = zeros(2*M - 1, N);

for j=2:length(B_x)-2
    for i=2:(N-1)
        % B_x = del(Az)/del(y)
        B_x(j, i) = (Az(j+1, i) - Az(j-1, i))/(2*(dy*0.01));
        B_y(j, i) = -((Az(j, i+1) - Az(j, i-1))/(2*(dx*0.01)));
    end
end

step = 1/dx;

%q = quiver(X(1:step:end, 1:step:end), Y(1:step:end, 1:step:end), B_x(1:step:end, 1:step:end), B_y(1:step:end, 1:step:end));
%q.AutoScaleFactor = 1.4;

mi = zeros(2*M-1, N);

mi_0 = 4*pi*10^-7;
mi_ferro = 2500*mi_0;
mi_ar = 1*mi_0;
mi_b = mi_ar;

for j=2:2*M-2
    for i=2:(N-1)
        if ((i-1)*dx < 4)
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

H_x = mi.*B_x;
H_y = mi.*B_y;

%q_h = quiver(X(1:step:end, 1:step:end), Y(1:step:end, 1:step:end), H_x(1:step:end, 1:step:end), H_y(1:step:end, 1:step:end));
%q_h.AutoScaleFactor = 1.4;


%{ 
q = quiver(X(1:step:end, 1:step:end), Y(1:step:end, 1:step:end), curl_x(1:step:end, 1:step:end), curl_y(1:step:end, 1:step:end));
q.AutoScaleFactor = 1.4;
%}

%}

% Item d)

f_x = zeros(length(Az));
f_y = zeros(length(Az));

f_x = (B_x(5/dx, :)).^2 + (B_y(5/dx, :)).^2;
f_y = 2*B_x(5/dx, :).*B_y(5/dx, :);

a = 1;
b = 2*M - 1;
n = length(Az);
h = dy;
s = f_x(a) + f_x(b);

for i = 1:2:n-1
    s = s + 4*f_x(a+i*h);
end

for i = 2:2:n-2
    s = s + 2*f_x(a+i*h);
end

I_x = (h*0.01/3) * s

s_y = f_y(a) + f_y(b);

for i = 1:2:n-1
    s_y = s_y + 4*f_y(a+i*h);
end

for i = 2:2:n-2
    s_y = s_y + 2*f_y(a+i*h);
end

I_y = (h*0.01/3) * s
