function plotGraph(t, estados, f_hist, metodo, h, item, caso)
%plotGraph: Imprime a resposta de theta e theta_dot em um gráfico e x_dot e
% x_dot2 em outro gráfico, todos em uma mesma figura.
%Input:
%   t       = vetor do tempo.
%   estados = vetor com os valores calculados em todas as iterações do
%             método utilizado.
%   f_hist  = vetor que contém os valores de x_dot2 calculados em todas as
%             iterações do método utilizado.
%   metodo  = nome do método utilizado.
%   h       = valor do passo considerado.
%   item    = número do item correspondente.
%   caso    = caso correspondente se o item for o 2, caso contrário não é
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
    title(['Método: ', metodo, ', h = ', num2str(h)], 'FontWeight', 'Normal', 'FontSize', 17);
else
    title({['Caso ', num2str(caso)], ['Método: ', metodo, ', h = ', num2str(h)]}, 'FontWeight', 'Normal', 'FontWeight', 'Normal', 'FontSize', 17);
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

