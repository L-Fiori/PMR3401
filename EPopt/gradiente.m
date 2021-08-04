function dVdxy = gradiente(mesh, V)

x = mesh.coor(mesh.con(1, :), 1);
y = mesh.coor(mesh.con(1, :), 2);

B = calc_B(x, y);

Vn = V(mesh.con');
dVdxy = (B*Vn)';


end

