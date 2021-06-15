clc;
clear all;
close all;

% Constantes
m_tot = 1939; % (kg)
L     = 2.95; % (m)
I     = 1;    % (kg * m^2)
omega = 10;   % (rad/s)
mi    = 0.42;
beta  = 0.02;
g     = 9.8;

% Item 1
% Condições iniciais

x0         = 0;
x_dot0     = 0;
theta0     = deg2rad(10);
theta_dot0 = 0;

m     = 0.6*m_tot;
m_1   = 0.4*m_tot;
F     = beta*m*g;
F_1   = mi*m_1*g;

t_0 = 0;
t_f = 20;
h = 1;

t = t_0:h:t_f;

estados0 = [theta_dot0; x_dot0; theta0; x0];

% Método de Euler

[estados, f_hist] = Euler(estados0, h, t, m_tot, m_1, L, F, F_1, I, omega);
plotGraph(t, estados, f_hist, 'Euler', h, 1, 0)

% RK2

[estados, f_hist] = RK2(estados0, h, t, m_tot, m_1, L, F, F_1, I, omega);
plotGraph(t, estados, f_hist, 'RK2', h, 1, 0)

% RK4

[estados, f_hist] = RK4(estados0, h, t, m_tot, m_1, L, F, F_1, I, omega);
plotGraph(t, estados, f_hist, 'RK4', h, 1, 0)

% Item 2

h = 1;
t = t_0:h:t_f;

% Caso 1
m   = 0.8*m_tot; % Motor dianteiro
m_1 = 0.2*m_tot; % Motor dianteiro
F   = beta*m*g;  % Tração traseira
F_1 = mi*m*g;    % Tração traseira

[estados, f_hist] = RK4(estados0, h, t, m_tot, m_1, L, F, F_1, I, omega);
plotGraph(t, estados, f_hist, 'RK4', h, 2, 1)

% Caso 2

m   = 0.8*m_tot;   % Motor dianteiro
m_1 = 0.2*m_tot;   % Motor dianteiro
F   = -mi*m*g;     % Tração dianteira
F_1 = -beta*m_1*g; % Tração dianteira

[estados, f_hist] = RK4(estados0, h, t, m_tot, m_1, L, F, F_1, I, omega);
plotGraph(t, estados, f_hist, 'RK4', h, 2, 2)

% Caso 3

m   = 0.2*m_tot;   % Motor traseiro
m_1 = 0.8*m_tot;   % Motor traseiro
F   = -mi*m*g;     % Tração dianteira
F_1 = -beta*m_1*g; % Tração dianteira

[estados, f_hist] = RK4(estados0, h, t, m_tot, m_1, L, F, F_1, I, omega);
plotGraph(t, estados, f_hist, 'RK4', h, 2, 3)

% Caso 4

m   = 0.2*m_tot; % Motor traseiro
m_1 = 0.8*m_tot; % Motor traseiro
F   = beta*m*g;  % Tração traseira
F_1 = mi*m*g;     % Tração traseira

[estados, f_hist] = RK4(estados0, h, t, m_tot, m_1, L, F, F_1, I, omega);
plotGraph(t, estados, f_hist, 'RK4', h, 2, 4)

% Caso 5a

m   = 0.2*m_tot; % Motor traseiro
m_1 = 0.8*m_tot; % Motor traseiro
F   = -mi*m*g;   % Tração 4 rodas
F_1 = mi*m_1*g;  % Tração 4 rodas

[estados, f_hist] = RK4(estados0, h, t, m_tot, m_1, L, F, F_1, I, omega);
plotGraph(t, estados, f_hist, 'RK4', h, 2, 5.1)

% Caso 5b

m   = 0.8*m_tot; % Motor dianteiro
m_1 = 0.2*m_tot; % Motor dianteiro
F   = -mi*m*g;   % Tração 4 rodas
F_1 = mi*m_1*g;  % Tração 4 rodas

[estados, f_hist] = RK4(estados0, h, t, m_tot, m_1, L, F, F_1, I, omega);
plotGraph(t, estados, f_hist, 'RK4', h, 2, 5.2)
