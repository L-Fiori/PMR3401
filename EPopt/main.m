clear; clc; close all;

Lx = 0.32;
Ly = 0.3;
nelx = 16;
nely = 15;

mesh = generate_mesh(Lx, Ly, nelx, nely);

% Plot da malha com os centroides
Nel = size(mesh.con, 1);
for e=1:Nel
	xe = mesh.coor(mesh.con(e,:), 1);
	ye = mesh.coor(mesh.con(e,:), 2);
	plot([xe; xe(1)], [ye; ye(1)], '-b');
	hold on;
end

plot(mesh.cen(:, 1), mesh.cen(:, 2), '*r');
axis equal;

Ngdl = size(mesh.coor, 1);

% Constantes de condutividade
sigma_a = 5.9e7;
sigma_b = 6.3e7;
sigma_c = 4.6e7;
sigma_h = 1e-8;

% Comeco atribuindo a condutividade do material B para toda a malha
sigma = sigma_b*ones(Nel, 1);

% Identificar os elementos de cada material
list = (1:Nel)';

id1a = ( mesh.cen(:, 2) <= 0.08 + mesh.cen(:, 1) <= 0.22 + mesh.cen(:, 1) >= 0.1 ) == 3;
id2a = ( mesh.cen(:, 2) >= 0.1 + mesh.cen(:, 2) <= 0.18 ) == 2;
id_a = [id1a id2a];

id1c = ( mesh.cen(:, 1) >= 0.1 + mesh.cen(:, 1) <= 0.22 + mesh.cen(:, 2) >= 0.26 ) == 3;
id2c = ( mesh.cen(:, 1) >= 0.06 + mesh.cen(:, 1) <= 0.26 + mesh.cen(:, 2) >= 0.22 + mesh.cen(:, 2) <= 0.26 ) == 4;
id3c = ( mesh.cen(:, 2) >= 0.08 + mesh.cen(:, 2) <= 0.12 ) == 2; 
id_c = [id1c id2c id3c];

% Aplicar as condutividades de cada material
sigma(list(id_a)) = sigma_a;
sigma(list(id_c)) = sigma_c;

% Identificar os elementos relativos ao buraco
id1h = ( mesh.cen(:, 1) >= 0.06 + mesh.cen(:, 1) <= 0.26 + mesh.cen(:, 2) >= 0.12 + mesh.cen(:, 2) <= 0.16 ) == 4;
id2h = ( mesh.cen(:, 1) >= 0.1 + mesh.cen(:, 1) <= 0.22 + mesh.cen(:, 2) >= 0.16 + mesh.cen(:, 2) <= 0.2 ) == 4;
id_h = [id1h id2h];

% Aplicar condutividade proxima de zero para os elementos do buraco
sigma(list(id_h)) = sigma_h;

% Plots de cada regiao da malha
plot(mesh.cen(:, 1), mesh.cen(:, 2), '*r'); hold on;

plot(mesh.cen(id1a, 1), mesh.cen(id1a, 2), 'ob'); hold on;
plot(mesh.cen(id2a, 1), mesh.cen(id2a, 2), 'ob'); hold on;

plot(mesh.cen(id1c, 1), mesh.cen(id1c, 2), 'sb'); hold on;
plot(mesh.cen(id2c, 1), mesh.cen(id2c, 2), 'sb'); hold on;
plot(mesh.cen(id3c, 1), mesh.cen(id3c, 2), 'sb'); hold on;

plot(mesh.cen(id1h, 1), mesh.cen(id1h, 2), 'og'); hold on;
plot(mesh.cen(id2h, 1), mesh.cen(id2h, 2), 'og'); hold on;

axis equal;

% Montagem da matriz global (usando for)
KG = zeros(Ngdl, Ngdl);
Ke = k_ele(mesh, 1);

for e=1:Nel
	id = mesh.con(e, :);
	KG(id, id) = KG(id, id) + sigma(e)*Ke;
end

% Montagem da matriz global (indexada)
Ke = k_ele(mesh, 1);
I = reshape(repmat(mesh.con, 1, 4)', Nel*16, 1);
J = mesh.con';
J = kron( J(:), ones(4, 1));
k_glob = repmat(Ke(:), Nel, 1).*kron(sigma, ones(16, 1));
KG = sparse(I, J, k_glob);

% Condicoes de contorno

list = (1:Ngdl)';
id1 = mesh.coor(:, 2) <= dy/2;
id2 = mesh.coor(:, 2) >= (Ly - dy/2);

V_cont = [list(id1) 0*ones(sum(id1), 1) ;
          list(id2) 250*ones(sum(id2), 1)];

[KGM, FM] = set_boundaryconditions(KG, V_cont);

V = KGM\FM;

% Pos-processamento



