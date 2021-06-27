function [B_x, B_y, H_x, H_y] = itemC(Az, B_x, B_y, dx, dy, N, mi_x, mi_y, isE) 

	for j=2:length(B_x)-1
	    for i=2:(N-1)
		% B_x = del(Az)/del(y)
		B_x(j, i) = (Az(j+1, i) - Az(j-1, i))/(2*dy*0.01);
		
		if (i-1)*deltax == 4 % Borda direita da armadura
		    By(j,i) = -(Az(j, i-2) - 4*Az(j, i-1) + 3*Az(j, i))/(2*dx*0.01); % Diferença regressiva de três pontos
		else
		    B_y(j, i) = -((Az(j, i+1) - Az(j, i-1))/(2*dx*0.01));
		end

	    end
	end

	if isE
		H_x = mi_y.*B_x;
	else
		H_x = mi_x.*B_x;
	end

	H_y = mi_y.*B_y;
end