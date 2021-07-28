function plot_struct(coorNormal, coords, con, str)

Nel = size(con, 1);

for e=1:Nel
    
    x1N = coorNormal(con(e, 1), 1);
	x2N = coorNormal(con(e, 2), 1);

	y1N = coorNormal(con(e, 1),2);
	y2N = coorNormal(con(e, 2),2);
    
    plot([x1N x2N], [y1N y2N], '--b'); hold on;

	x1 = coords(con(e, 1), 1);
	x2 = coords(con(e, 2), 1);

	y1 = coords(con(e, 1),2);
	y2 = coords(con(e, 2),2);

	plot([x1 x2], [y1 y2], str); hold on;

end
axis equal;
title('Primeiro modo de vibrar')

end