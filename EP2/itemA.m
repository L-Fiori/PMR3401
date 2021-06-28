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

iter = true; 
count = 0; 
lambda = 1.75;
E = 0.0001;

front_v = false;
front_h = false;

while iter
    highest_E = 0;
    for j = 1:((10/dy) + 1)
        for i = 2:((22/dx) + 1)
            front_v = false;
            front_h = false;

            Az_old = Az(j, i);
            
            % Caso o ponto pertenca ao contorno
            if((i-1)*dx == 0 || (i-1)*dx == 22 || (j-1)*dy == 10 || ((i-1)*dx == 20 && ((j-1)*dy >= 6 && (j-1) < 10) || ((i-1)*dx > 20 && (i-1)*dx < 22 && (j-1)*dy == 6)))
                Az(j, i) = 0;
                
            % Caso o ponto pertenca a bobina lado esquerdo
            elseif(((i-1)*dx > 14 && (i-1)*dx < 16) && ((j-1)*dy < 6))
                if(j == 1)
                    Az(j,i) = (2*Az(j+1, i) + Az(j, i+1) + Az(j, i-1) - ((dx*0.01)^2)*mi_b*(2*(10^6)*cos((pi*(j-1)/dy)*0.01/0.12) + 8*10^5)) / 4;
                else
                    Az(j,i) = (Az(j+1, i) + Az(j-1, i) + Az(j, i+1) + Az(j, i-1) - ((dx*0.01)^2)*mi_b*(2*(10^6)*cos((pi*(j-1)/dy)*0.01/0.12) + 8*10^5)) / 4;
                end
                
            % Caso o ponto pertenca a bobina lado direito
            elseif(((i-1)*dx > 20 && (i-1)*dx < 22) && ((j-1)*dy < 6))
                if(j == 1)
                    Az(j,i) = (2*Az(j+1, i) + Az(j, i+1) + Az(j, i-1) + ((dx*0.01)^2)*mi_b*(2*(10^6)*cos((pi*(j-1)/dy)*0.01/0.12) + 8*10^5)) / 4;
                else
                    Az(j,i) = (Az(j+1, i) + Az(j-1, i) + Az(j, i+1) + Az(j, i-1) + ((dx*0.01)^2)*mi_b*(2*(10^6)*cos((pi*(j-1)/dy)*0.01/0.12) + 8*10^5)) / 4;
                end
                
            elseif((i-1)*dx == 4)
                % Fronteira vertical ferro-ar
                front_v = true;
                mi_2 = mi_ar;
                mi_1 = mi_ferro;

            elseif((i-1)*dx == 5 && (j-1)*dy >= 6)
                % Fronteira vertical ar-ferro
                front_v = true;
                mi_2 = mi_ferro;
                mi_1 = mi_ar;
                
            elseif(((i-1)*dx > 5 && (i-1)*dx < 14) && (j-1)*dy == 6)
                % Fronteira horizontal ar-ferro
                front_h = true;
                mi_2 = mi_ferro;
                mi_1 = mi_ar;
                
            elseif(((i-1)*dx >= 14 && (i-1)*dx <= 16) && (j-1)*dy == 6)
                % Fronteira horizontal bobina-ferro
                front_h = true;
                mi_2 = mi_ferro;
                mi_1 = mi_b;
            
            elseif(((i-1)*dx == 14 && (j-1)*dy < 6))
                % Fronteira vertical ar-bobina
                front_v = true;
                mi_2 = mi_b;
                mi_1 = mi_ar;
                
            elseif(((i-1)*dx == 16 && (j-1)*dy < 6))
                % Fronteira vertical bobina-ferro
                front_v = true;
                mi_2 = mi_ferro;
                mi_1 = mi_b;
            
            elseif(((i-1)*dx == 20 && (j-1)*dy < 6))
                % Fronteira vertical ferro-bobina
                front_v = true;
                mi_2 = mi_b;
                mi_1 = mi_ferro;
            
            % Caso seja um ponto interno do dominio sem ser fronteira ou bobina (a maioria dos pontos)
            else
                if(j == 1)
                    Az(j, i) = (2*Az(j+1, i) + Az(j, i+1) + Az(j, i-1))/4;
                else
                    Az(j, i) = (Az(j+1, i) + Az(j-1, i) + Az(j, i+1) + Az(j, i-1))/4;
                end
            end
           
            if front_v
                if(j == 1)
                    Az(j, i) = (1/(4*((mi_2/mi_1) + 1)))*(((mi_2/mi_1) + 1)*2*Az(j+1, i) + 2*(Az(j, i+1) + (mi_2/mi_1)*Az(j, i-1)));
                else
                    Az(j, i) = (1/(4*(mi_2/mi_1+1)))*((mi_2/mi_1+1)*(Az(j+1, i) + Az(j-1, i)) + 2*(Az(j, i+1) + (mi_2/mi_1)*Az(j, i-1)));
                end
            end

            if front_h
                if(j == 1)
                    Az(j, i) = (1/(4*((mi_2/mi_1) + 1)))*(((mi_2/mi_1) + 1)*(Az(j, i+1) + Az(j, i-1)) + 2*(Az(j+1, i) + (mi_2/mi_1)*Az(j+1, i)));
                else
                    Az(j, i) = (1/(4*((mi_2/mi_1) + 1)))*(((mi_2/mi_1) + 1)*(Az(j, i+1) + Az(j, i-1)) + 2*(Az(j+1, i) + (mi_2/mi_1)*Az(j-1, i)));
                end
            end

            Az(j,i) = lambda*Az(j,i) + (1-lambda)*Az_old;
            % Armazena o maior valor do erro em modulo
            current_E = abs((Az(j,i) - Az_old)/Az(j, i));
            if (current_E > highest_E)
                highest_E = current_E;
            end 

        end
    end

    if highest_E < E
        iter = false;
    end

    count = count + 1;

end

Az

display('Iteracoes: ')
count

end