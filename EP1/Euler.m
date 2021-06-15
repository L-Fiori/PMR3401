function [estados, f_hist] = Euler(estados, h, t, m_tot, m_1, L, F, F_1, I, omega)
%Euler: Resolve uma EDO genérica através do método de Euler.
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

for n=1:length(t)-1
    derivadas = f(estados(1, n), estados(3, n), m_tot, m_1, L, F, F_1, I, omega);
    K1(1:2) = derivadas; % theta_dot2 (i) e x_dot2 (i)
    f_hist(n) = K1(2);
    K1(3:4) = estados(1:2, n); % theta_dot (i) e x_dot (i)
    estados(:, n+1) = estados(:, n) + (h * K1);
end

    derivadas = f(estados(1, length(estados)), estados(3, length(estados)), m_tot, m_1, L, F, F_1, I, omega);
    f_hist(length(f_hist) + 1) = derivadas(2);

end


