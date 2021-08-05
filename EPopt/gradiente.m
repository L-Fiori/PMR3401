function dVdxy = gradiente(mesh, V)

Nel = size(mesh.con, 1);
dVdxy = zeros(Nel, 2);

for e=1:Nel
    x = mesh.coor(mesh.con(e, :), 1);
    y = mesh.coor(mesh.con(e, :), 2);

    B = calc_B(x, y);

    Vn = V(mesh.con(e, :)');
    
    dVdxy(e, :) = (B*Vn)';
end

