function plotGraph(t, estados, f_hist, metodo, h, item, caso)
%plotGraph: Imprime a resposta de theta e theta_dot em um gr�fico e x_dot e
% x_dot2 em outro gr�fico, todos em uma mesma figura.
%Input:
%   t       = vetor do tempo.
%   estados = vetor com os valores calculados em todas as itera��es do
%             m�todo utilizado.
%   f_hist  = vetor que cont�m os valores de x_dot2 calculados em todas as
%             itera��es do m�todo utilizado.
%   metodo  = nome do m�todo utilizado.
%   h       = valor do passo considerado.
%   item    = n�mero do item correspondente.
%   caso    = caso correspondente se o item for o 2, caso contr�rio n�o �
%             utilizado.

figure

subplot(2, 1, 1);
yyaxis left
plot(t, estados(1, :));
ylabel('$\dot{\theta} (rad/s)$')
yyaxis right
plot(t, estados(3, :));
ylabel('$\theta (rad)$');
xlabel('t (s)')
if item == 1
    title(['M�todo: ', metodo, ', h = ', num2str(h)], 'FontWeight', 'Normal', 'FontSize', 17);
else
    title({['Caso ', num2str(caso)], ['M�todo: ', metodo, ', h = ', num2str(h)]}, 'FontWeight', 'Normal', 'FontWeight', 'Normal', 'FontSize', 17);
end
grid on;

subplot(2, 1, 2);
yyaxis left
plot(t, estados(2, :));
ylabel('$\dot{x} (m/s)$')
yyaxis right
plot(t, f_hist);
ylabel('$\ddot{x} (m/s^2)$')
xlabel('t (s)')
grid on;


end

