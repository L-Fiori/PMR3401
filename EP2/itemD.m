function [F_x, F_y] = itemD(B_x, B_y, dx, M, n, mi_0) 

	f_x = (B_x((4/dx)+1, :)).^2 - (B_y((4/dx)+1, :)).^2
	f_y = 2*B_x((4/dx)+1, :).*B_y((4/dx)+1, :)

	B_x((4/dx)+1, :)
	B_y((4/dx)+1, :)

	%B_x(5/dx, :);
	%f_x;
	%f_y;

	a = 1;
	b = 2*M - 1;
	h = dx;
	s_x = f_x(a) + f_x(b);

	for i = 1:2:n-1
	s_x = s_x + 4*f_x(a+i);
	end

	for i = 2:2:n-2
	s_x = s_x + 2*f_x(a+i);
	end

	F_x = (1/(2*mi_0))*(h*0.01/3) * s_x

	s_y = f_y(a) + f_y(b);

	for i = 1:2:n-1
	s_y = s_y + 4*f_y(a+i);
	end

	for i = 2:2:n-2
	s_y = s_y + 2*f_y(a+i);
	end

	F_y = (1/(2*mi_0))*(h*0.01/3) * s_y

end