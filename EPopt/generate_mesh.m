function mesh = generate_mesh(Lx, Ly, nelx, nely)

nodx = nelx + 1;
nody = nely + 1;

x = linspace(0, Lx, nodx)';
y = linspace(0, Ly, nody)';
x = repmat(x, nody, 1);
y = kron(y, ones(nodx, 1));
mesh.coor = [x y];

Nel = nelx*nely;
mesh.con = zeros(Nel, 3);
cont = 0;
for j=1:nely
	for i=1:nelx
		cont = cont + 1;
		idx = i + (j-1)*nodx;
		idy = i + j*nodx;
		mesh.con((2*cont-1):(2*cont), :) = [idx idx+1 idy; idx+1 idy+1 idy];
	end
end

% Criacao dos centroides
mesh.cen = zeros(2*Nel, 2);
for e=1:2*Nel
    mesh.cen(e, :) = [mean(mesh.coor(mesh.con(e,:), 1)) mean(mesh.coor(mesh.con(e,:), 2))];
end
    
%dx = Lx/nelx;
%dy = Ly/nely;

%x = linspace(dx/2, Lx-dx/2, nelx)';
%y = linspace(dy/2, Ly-dy/2, nely)';
%x = repmat(x, nely, 1);

end