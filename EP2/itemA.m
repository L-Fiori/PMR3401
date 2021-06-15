function [Az] = ItemA(Az, dx, dy)
%Input:
%   Az =
%   dx =
%   dy = 
%Output:
%   Az = 

mi_0 = 4*pi*10^-7;
mi_ferro = 2500*mi_0;
mi_ar = 1*mi_0;
mi_b = mi_ar;

iter = true; % so para de iterar quando iter for false
count = 0; % Contador de iteracoes
lambda = 1.75;
E = 0.0001;

while iter == true
    highest_E = 0;
    for j = 1:(10/dy)
        for i = 2:(22/dx)
            Az_old = Az(j, i);
            
            % Caso o ponto pertenca ao contorno
            if(i*dx == 0 || i*dx == 22 || j*dy == 10 || (i*dx == 20 &&(j*dy >= 6 && j < 10) || (i*dx > 20 && i*dx < 22 && j*dy == 6)))
                continue
                
            % Caso o ponto pertenca a bobina
            elseif(((i*dx > 14 && i*dx < 16) || (i*dx > 20 && i*dx < 22)) && (j*dy < 6))
                if(j == 1)
                    Az(j,i) = (2*Az(j+1, i) + Az(j, i+1) + Az(j, i-1) + (dx^2)*mi_b*(2*(10^6)*cos((pi*i/dy)/0.12) + 8*10^5)) / 4;
                else
                    Az(j,i) = (Az(j+1, i) + Az(j-1, i) + Az(j, i+1) + Az(j, i-1) + (dx^2)*mi_b*(2*(10^6)*cos((pi*i/dy)/0.12) + 8*10^5)) / 4;
                end

            elseif(i*dx == 4)
                % Fronteira vertical ferro-ar
                if(j == 1)
                    Az(j, i) = (1/(4*((mi_ar/mi_ferro) + 1)))*(((mi_ar/mi_ferro) + 1)*2*Az(j+1, i) + 2*(Az(j, i+1) + (mi_ar/mi_ferro)*Az(j, i-1)));
                else
                    Az(j, i) = (1/(4*((mi_ar/mi_ferro) + 1)))*(((mi_ar/mi_ferro) + 1)*(Az(j+1, i) + Az(j-1, i)) + 2*(Az(j, i+1) + (mi_ar/mi_ferro)*Az(j, i-1)));
                end
 
            elseif(i*dx == 5 && j*dy >= 6)
                % Fronteira vertical ar-ferro
                if(j*dy == 6)
                    Az(j, i) = (1/(4*((mi_ferro/mi_ar) + 1))) * ((mi_ferro/mi_ar)*(Az(j-1, i) + Az(j, i-1)) + Az(j+1, i) + Az(j, i+1));
                else 
                    Az(j, i) = (1/(4*((mi_ferro/mi_ar) + 1)))*(((mi_ferro/mi_ar) + 1)*(Az(j+1, i) + Az(j-1, i)) + 2*(Az(j, i+1) + (mi_ferro/mi_ar)*Az(j, i-1)));
                end

            elseif((i*dx > 5 && i*dx < 14) && j*dy == 6)
                % Fronteira horizontal ar-ferro
                Az(j, i) = (1/(4*((mi_ferro/mi_ar) + 1))) * (((mi_ferro/mi_ar) + 1)*(Az(j, i+1) + Az(j, i-1)) + 2*(Az(j+1, i) + (mi_ferro/mi_ar)*Az(j-1, i)));
            
            elseif((i*dx >= 14 && i*dx <= 16) && j == 6)
                % Fronteira horizontal bobina-ferro
                if(j == 1)
                    Az(j, i) = (1/(4*((mi_ferro/mi_b) + 1))) * (((mi_ferro/mi_b) + 1)*(Az(j, i+1) + Az(j, i-1)) + 2*(Az(j+1, i) + (mi_ferro/mi_b)*Az(j+1, i)));
                else
                    Az(j, i) = (1/(4*((mi_ferro/mi_b) + 1))) * (((mi_ferro/mi_b) + 1)*(Az(j, i+1) + Az(j, i-1)) + 2*(Az(j+1, i) + (mi_ferro/mi_b)*Az(j-1, i)));
                end
            
            elseif((i*dx == 14 && j*dy < 6))
                % Fronteira vertical ar-bobina
                if (j == 1)
                    Az(j, i) = (1/(4*((mi_b/mi_ar) + 1))) * ((mi_b/mi_ar)*(Az(j+1, i) + Az(j, i-1)) + Az(j+1, i) + Az(j, i+1));
                else 
                    Az(j, i) = (1/(4*((mi_b/mi_ar) + 1)))*(((mi_b/mi_ar) + 1)*(Az(j+1, i) + Az(j-1, i)) + 2*(Az(j, i+1) + (mi_b/mi_ar)*Az(j, i-1)));
                end

            elseif((i*dx == 16 && j*dy < 6))
                % Fronteira vertical bobina-ferro
                if (j == 1)
                    Az(j, i) = (1/(4*((mi_ferro/mi_b) + 1))) * ((mi_ferro/mi_b)*(Az(j+1, i) + Az(j, i-1)) + Az(j+1, i) + Az(j, i+1));
                else 
                    Az(j, i) = (1/(4*((mi_ferro/mi_b) + 1)))*(((mi_ferro/mi_b) + 1)*(Az(j+1, i) + Az(j-1, i)) + 2*(Az(j, i+1) + (mi_ferro/mi_b)*Az(j, i-1)));
                end
            
            elseif((i*dx == 20 && j*dy < 6))
                % Fronteira vertical ferro-bobina
                if (j == 1)
                    Az(j, i) = (1/(2*((mi_b/mi_ferro) + 1))) * ((mi_b/mi_ferro)*(Az(j+1, i) + Az(j, i-1)) + Az(j+1, i) + Az(j, i+1));
                else 
                    Az(j, i) = (1/4*((mi_b/mi_ferro) + 1))*(((mi_b/mi_ferro) + 1)*(Az(j+1, i) + Az(j-1, i)) + 2*(Az(j, i+1) + (mi_b/mi_ferro)*Az(j, i-1)));
                end
            
                % Caso seja um ponto interno do dominio sem ser fronteira ou bobina (a maioria dos pontos)
            else
                if(j == 1)
                    Az(j, i) = (2*Az(j+1, i) + Az(j, i+1) + Az(j, i-1))/4;
                else
                    Az(j, i) = (Az(j+1, i) + Az(j-1, i) + Az(j, i+1) + Az(j, i-1))/4;
                end
            end
            
            Az(j,i) = lambda*Az(j,i) + (1-lambda)*Az_old;
            % Armazena o maior valor do erro em modulo
            current_E = abs(Az(j,i) - Az_old);
            if (current_E > highest_E)
                j
                i
                highest_E = current_E
            end 

        end
    end

    if highest_E < E
        iter = false;
    end

    count = count + 1;

end

display('Iteracoes: ')
count

end