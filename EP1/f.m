function [derivadas] = f(theta_dot, theta, m_tot, m_1, L, F, F_1, I, omega)
%f: Calcula as derivadas theta_dot2 e x_dot2, oriundas de manipula��o 
%   alg�brica das equa��es provenientes da modelagem do problema.
%Input:
%   theta_dot = primeira derivada temporal de theta.
%   theta     = �ngulo theta considerado na modelagem do problema.
%   m_tot, m_1, L, F, F_1, I s�o constantes.
%Output:
%   derivadas = vetor que cont�m o resultado do c�lculo das derivadas
%               theta_dot2 e x_dot2 nessa ordem.

theta_dot2 = (sin(theta)/(L*((sin(theta)^2)*m_1 - m_tot)))*(F_1*cos(theta) - F - m_1*L*(theta_dot^2)*cos(theta)+ ((2*I*omega*theta_dot*m_tot)/(m_1*L*sin(theta))));
x_dot2 = (F_1*cos(theta) - F - m_1*L*(theta_dot^2)*cos(theta) + 2*I*omega*theta_dot*sin(theta)/L)/(m_tot - m_1*(sin(theta)^2));

derivadas = [theta_dot2; x_dot2];

end

