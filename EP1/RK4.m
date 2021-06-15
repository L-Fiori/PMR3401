function [estados, f_hist] = RK4(estados, h, t, m_tot, m_1, L, F, F_1, I, omega)
%RK4: Resolve uma EDO genérica através do método de Runge-Kutta de
%     quarta ordem.
%Input:
%   estados = vetor que contém os valores calculados na iteração 
%             anterior (ou os valores iniciais para o caso de ser a primeira iteração),
%             a saber: theta_dot, x_dot, theta e x nessa ordem.
%   h       = valor do passo considerado.
%   t       = vetor do tempo.
%   m_tot, m_1, L, F, F_1, I são constantes passadas como entrada para que
%   seja realizado o cálculo das derivadas de segunda ordem.
%Output:
%   estados = vetor com os valores calculados em todas as iterações do
%             método.
%   f_hist  = vetor que contém os valores de x_dot2 calculados em todas as
%             iterações do método.

K1 = [0; 0; 0; 0];
K2 = [0; 0; 0; 0];
K3 = [0; 0; 0; 0];
K4 = [0; 0; 0; 0];

for n=1:length(t)-1
    derivadas = f(estados(1, n), estados(3, n), m_tot, m_1, L, F, F_1, I, omega);
    K1(1:2) = derivadas; % theta_dot2 (i) e x_dot2 (i)
    f_hist(n) = K1(2);
    K1(3:4) = estados(1:2, n); % theta_dot (i) e x_dot (i)
    derivadas = f(estados(1,n) + (h/2)*K1(1), estados(3, n) + (h/2)*K1(3), m_tot, m_1, L, F, F_1, I, omega);
    K2(1:2) = derivadas;
    K2(3:4) = estados(1:2, n) + (h/2)*K1(3:4);
    derivadas = f(estados(1,n) + (h/2)*K2(1), estados(3, n) + (h/2)*K2(3), m_tot, m_1, L, F, F_1, I, omega);
    K3(1:2) = derivadas;
    K3(3:4) = estados(1:2, n) + (h/2)*K2(3:4);
    derivadas = f(estados(1,n) + h*K3(1), estados(3, n) + h*K3(3), m_tot, m_1, L, F, F_1, I, omega);
    K4(1:2) = derivadas;
    K4(3:4) = estados(1:2, n) + h*K3(3:4);
    
    K = K1 + 2*K2 + 2*K3 + K4;
    estados(:, n+1) = estados(:, n) + (h/6 * K);
end

    derivadas = f(estados(1, length(estados)), estados(3, length(estados)), m_tot, m_1, L, F, F_1, I, omega);
    f_hist(length(f_hist) + 1) = derivadas(2);

end

