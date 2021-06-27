function [Az_new] = itemE2(Az, dx, dt, t_range, divs, M, N)

sigma = 4*10^6;
t = dt;
tamanho = t_range/divs;

mi_0 = 4*pi*10^-7;
mi_ferro_x = 1200*mi_0;
mi_ferro_y = 2500*mi_0;
mi_ar = 1*mi_0;
mi_b = mi_ar;

Az_new = zeros(M, N, tamanho);
count = 1;

for k=1:t_range
    
    % y = 0 -> -10
    for j=1:(M-1) % Linhas

        % x = 0 -> 22
        for i=2:(N-1) % Colunas

            % ---- Contorno Externo ----
            % Já resolvido Az = 0

            % ---- Fronteiras Internas ----
            % Condicao de contorno de x = 4, 0 > y > -10 (Armadura/Ar em y)
            if ((i-1)*dx == 4)
                mi_1x = mi_ferro_x;
                mi_1y = mi_ferro_y;
                mi_2 = mi_ar;
                if (j == 1)
                    Az(j,i,1) = 1/(4+2*mi_2*(1/mi_1y+1/mi_1x)) * ((mi_2/mi_1x+1)*(2*Az(j+1,i,1)) + 2*Az(j,i+1,1) + 2*(mi_2/mi_1y)*Az(j,i-1,1));
                else
                    Az(j,i,1) = 1/(4+2*mi_2*(1/mi_1y+1/mi_1x)) * ((mi_2/mi_1x+1)*(Az(j+1,i,1)+Az(j-1,i,1)) + 2*Az(j,i+1,1) + 2*(mi_2/mi_1y)*Az(j,i-1,1));
                end
                Az(j,i,2) = Az(j,i,1);

            % Condicao de contorno de x = 5, -6 > y > -10 (Ar/Ferro em y)
            elseif ((i-1)*dx == 5 && (j-1)*dx >= 6)
                mi_1 = mi_ar;
                mi_2x = mi_ferro_x;
                mi_2y = mi_ferro_y;
                if ((j-1)*dx == 6) % Quina
                    % Forma calculada: ???
                    Az(j,i,1) = 1/(2*mi_2y*(2/mi_1+1/mi_2y+1/mi_2x)) * ( mi_2y*(1/mi_2x+1/mi_1)*(Az(j+1,i,1)+Az(j-1,i,1)) + 2*Az(j,i+1,1) + 2*(mi_2y/mi_1)*Az(j,i-1,1) );
                else
                    Az(j,i,1) = 1/(2*mi_2y*(2/mi_1+1/mi_2y+1/mi_2x)) * ( mi_2y*(1/mi_2x+1/mi_1)*(Az(j+1,i,1)+Az(j-1,i,1)) + 2*Az(j,i+1,1) + 2*(mi_2y/mi_1)*Az(j,i-1,1) );
                end
                Az(j,i,2) = Az(j,i,1);

            % Condicao de contorno de 5 < x < 16, y = -6 (Ar/Ferro em x)
            elseif ((j-1)*dx == 6 && (i-1)*dx > 5 && (i-1)*dx < 16)
                mi_1 = mi_ar;
                mi_2x = mi_ferro_x;
                mi_2y = mi_ferro_y;
                Az(j,i,1) = 1/(2*mi_2y*(2/mi_1+1/mi_2y+1/mi_2x)) * ( mi_2y*(1/mi_2x+1/mi_1)*(Az(j,i+1,1)+Az(j,i-1,1)) + 2*Az(j+1,i,1) + 2*(mi_2y/mi_1)*Az(j-1,i,1) );
                Az(j,i,2) = Az(j,i,1);
                
            % Condicao de contorno de x = 14, 0 > y > -6 (Ar/Bobina em y)
            % (Igual ponto interno)

            % Condicao de contorno de x = 16, 0 > y > -6 (Bobina/Ferro em y)
            elseif ((i-1)*dx == 16 && (j-1)*dx < 6)
                mi_1 = mi_b;
                mi_2x = mi_ferro_x;
                mi_2y = mi_ferro_y;
                if (j == 1)
                    Az(j,i,1) = 1/(2*mi_2y*(2/mi_1+1/mi_2y+1/mi_2x)) * ( mi_2y*(1/mi_2x+1/mi_1)*(2*Az(j+1,i,1)) + 2*Az(j,i+1,1) + 2*(mi_2y/mi_1)*Az(j,i-1,1) );
                else
                    Az(j,i,1) = 1/(2*mi_2y*(2/mi_1+1/mi_2y+1/mi_2x)) * ( mi_2y*(1/mi_2x+1/mi_1)*(Az(j+1,i,1)+Az(j-1,i,1)) + 2*Az(j,i+1,1) + 2*(mi_2y/mi_1)*Az(j,i-1,1) );
                end
                Az(j,i,2) = Az(j,i,1);

            % Condicao de contorno de x = 20, 0 > y > -6 (Ferro/Bobina em y)
            elseif ((i-1)*dx == 20 && (j-1)*dx < 6)
                mi_1x = mi_ferro_x;
                mi_1y = mi_ferro_y;
                mi_2 = mi_b;
                if (j == 1)
                    Az(j,i,1) = 1/(4+2*mi_2*(1/mi_1y+1/mi_1x)) * ( (mi_2/mi_1x+1)*(2*Az(j+1,i,1)) + 2*Az(j,i+1,1) + 2*(mi_2/mi_1y)*Az(j,i-1,1) );
                else
                    Az(j,i,1) = 1/(4+2*mi_2*(1/mi_1y+1/mi_1x)) * ( (mi_2/mi_1x+1)*(Az(j+1,i,1)+Az(j-1,i,1)) + 2*Az(j,i+1,1) + 2*(mi_2/mi_1y)*Az(j,i-1,1) );
                end
                Az(j,i,2) = Az(j,i,1);

            % --- Pontos Internos ---
            % Região da bobina
            elseif ((i-1)*dx < 20 || (j-1)*dx < 6) 
                % Dentro da bobina esquerda
                if ((i-1)*dx > 14 && (i-1)*dx < 16 && (j-1)*dx < 6)
                    if (j == 1)
                        Az(j,i,2) = dt/(mi_b*sigma*(dx*10^-2)^2) * ( Az(j,i,1)*(-4 + mi_b*sigma*(dx*10^-2)^2/dt) - mi_b*(dx*10^-2)^2*cos(60*t)*(2*10^6*cos(pi*(j-1)*dx*10^-2/0.12)+8*10^5)+Az(j,i+1,1)+Az(j,i-1,1)+2*Az(j+1,i,1) );
                    else
                        Az(j,i,2) = dt/(mi_b*sigma*(dx*10^-2)^2) * (Az(j,i,1)*(-4 + mi_b*sigma*(dx*10^-2)^2/dt) - mi_b*(dx*10^-2)^2*cos(60*t)*(2*10^6*cos(pi*(j-1)*dx*10^-2/0.12)+8*10^5) + Az(j,i+1,1) + Az(j,i-1,1) + Az(j+1,i,1) + Az(j-1,i,1) );
                    end
                % Dentro da bobina direita
                elseif ((i-1)*dx > 20 && (i-1)*dx < 22 && (j-1)*dx < 6)
                    if (j == 1)
                        Az(j,i,2) = dt/(mi_b*sigma*(dx*10^-2)^2) * ( Az(j,i,1)*(-4 + mi_b*sigma*(dx*10^-2)^2/dt) + mi_b*(dx*10^-2)^2*cos(60*t)*...
                            (2*10^6*cos(pi*(j-1)*dx*10^-2/0.12)+8*10^5) + Az(j,i+1,1)+Az(j,i-1,1)+2*Az(j+1,i,1) );
                    else
                        Az(j,i,2) = dt/(mi_b*sigma*(dx*10^-2)^2) * ( Az(j,i,1)*(-4 + mi_b*sigma*(dx*10^-2)^2/dt) + mi_b*(dx*10^-2)^2*cos(60*t)*...
                            (2*10^6*cos(pi*(j-1)*dx*10^-2/0.12)+8*10^5)+...
                            Az(j,i+1,1)+Az(j,i-1,1)+...
                            Az(j+1,i,1)+Az(j-1,i,1) );
                    end
                % Dentro da armadura ou do ferro
                elseif ((i-1)*dx < 4 || ((i-1)*dx > 5 && (j-1)*dx > 6) || ((i-1)*dx > 16 && (i-1)*dx < 20 && (j-1)*dx <= 6) || ((i-1)*dx == 16 && (j-1)*dx == 6) )
                    if (j == 1)
                        Az(j,i,1) = 1/(1/mi_ferro_x+1/mi_ferro_y) * ( (1/mi_ferro_y)*(Az(j,i+1,1)+Az(j,i-1,1)) + (1/mi_ferro_x)*(2*Az(j+1,i,1)) ) / 2;
                    else
                        Az(j,i,1) = 1/(1/mi_ferro_x+1/mi_ferro_y) * ( (1/mi_ferro_y)*(Az(j,i+1,1)+Az(j,i-1,1)) + (1/mi_ferro_x)*(Az(j+1,i,1)+Az(j-1,i,1)) ) / 2;
                    end
                    Az(j,i,2) = Az(j,i,1);
                % Demais pontos
                else
                    if (j == 1)
                        Az(j,i,1) = ( Az(j,i+1,1) + Az(j,i-1,1) + 2*Az(j+1,i,1) )/4;
                    else
                        Az(j,i,1) = ( Az(j,i+1,1) + Az(j,i-1,1) + Az(j-1,i,1) + Az(j+1,i,1) )/4;
                    end
                    Az(j,i,2) = Az(j,i,1);
                end
            end


            
        end
    end
    
    Az(:,:,1) = Az(:,:,2);
    
    if (mod(k, divs) == 0)
        Az_new(:,:,count) = Az(:,:,1);
        count = count+1;
    end
    
    t = t+dt;
    
end


end