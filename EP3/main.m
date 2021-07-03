A1 = 1.8*0.9;
A2 = 0.9*0.9;
D3 = 0.050;

% Passo 1 - Definir a geometria (coordenadas e conectividades).

if (dx == 1)
	coords = [ 5  5;
       		  18  5;
		  65  5;
		  36 23;
		  linspace(0, dx, 65)' 5;
		  36 linspace(0, dx, 23)'];

	Nod_vh = size(linspace(0, dx, 65));
	Nod_vv = size(linspace(0, dx, 23));
else
	coords = [ 5  5;
       		  18  5;
		  65  5;
		  36 23;
		  linspace(0, dx, 65)' 5;
		  65 5
		  36 linspace(0, dx, 23)'
		  36 23];

	Nod_vh = size(linspace(0, dx, 65)) + 1;
	Nod_vv = size(linspace(0, dx, 23)) + 1;
end

	con = [ 1 4;
		2 4;
		3 4;
		(5:(5 + Nod_vh - 2))' (6:(6 + Nod_vh - 2))';
		((5 + Nod_vh):(5 + Nod_vh + Nod_vv - 2))' ((6 + Nod_vh):(6 + Nod_vh + Nod_vv - 2))'];

plot_struct(coords, con, '-b');

% Passo 2 - Atribuir propriedades aos elementos.

Nod = size(coords, 1);
Nel = size(con, 1);
Nel_vh = Nod_vh - 1;
Nel_vv = Nod_vv - 1;
Ngdl = 2*Nod; % Talvez modifique porque teoricamente a viga vertical é um elemento de portico e elemento de portico tem 3 graus de liberdade.


data.E = 210e9*ones(Nel, 1); % Vetor de modulos de Young dos elementos.
data.rho = 7600*ones(Nel, 1); % Vetor de massa especifica dos elementos.
data.A = ones(Nel, 1); % Vetor de areas dos elementos.

data.A(1:3) = pi*(D3/2)^2;
data.A(4:(4 + Nel_vh - 1)) = A2;
data.A((4 + Nel_vh):(4 + Nel_vh + Nel_vv - 1)) = A1;

data.L = zeros(Nel, 1); % Vetor de comprimento dos elementos.
data.Q = zeros(Nel, 1); % Vetor de rotacoes dos elementos.

for e=1:Nel

	x1 = coords(con(e, 1), 1);
	x2 = coords(con(e, 2), 1);

	y1 = coords(con(e, 1),2);
	y2 = coords(con(e, 2),2);

	data.L(e) = sqrt( (x2 - x1)^2 + (y2 - y1)^2 );
	data.Q(e) = atan2((y2 - y1), (x2 - x1))*(180/pi);

end

% Matrizes globais de rigidez e massa.

Kg = zeros(Ngdl, Ngdl);
Mg = zeros(Ngdl, Ngdl);

% Construcao das matrizes para os elementos de trelica

k0 = [1  0 -1 0;
      0  0  0 0;
      -1 0  1 0;
      0  0  0 0];

m0 = [2 0 1 0;
      0 0 0 0;
      1 0 2 0;
      0 0 0 0];

for e=1:3

	% Passo 3 - Construir a matriz de cada elemento.
	ke = (data.E(e)*data.A(e)/data.L(e))*k0;
	me = (data.rho(e)*data.A(e)*data.L(e)/6)*m0;

	% Passo 4 - Rotacionar cada elemento para o sistema global.
	c = cosd(data.Q(e)); s = sind(data.Q(e));
	T = [c  s  0 0;
	     -s c  0 0;
	     0  0  c s;
	     0  0 -s c];
	ke = T'*ke*T;
	me = T'*me*T;

	% Passo 5 - Alocar a informacao de cada matriz de elemento na matriz global.
	nod1 = con(e, 1);
	nod2 = con(e, 2);

	Kg(2*nod1-1:2*nod1, 2*nod1-1:2*nod1) = Kg(2*nod1-1:2*nod1, 2*nod1-1:2*nod1) + ke(1:2, 1:2); 
	Kg(2*nod1-1:2*nod1, 2*nod2-1:2*nod2) = Kg(2*nod1-1:2*nod1, 2*nod2-1:2*nod2) + ke(1:2, 3:4); 
	Kg(2*nod2-1:2*nod2, 2*nod1-1:2*nod1) = Kg(2*nod2-1:2*nod2, 2*nod1-1:2*nod1) + ke(3:4, 1:2); 
	Kg(2*nod2-1:2*nod2, 2*nod2-1:2*nod2) = Kg(2*nod2-1:2*nod2, 2*nod2-1:2*nod2) + ke(3:4, 3:4); 
	
	Mg(2*nod1-1:2*nod1, 2*nod1-1:2*nod1) = Mg(2*nod1-1:2*nod1, 2*nod1-1:2*nod1) + me(1:2, 1:2); 
	Mg(2*nod1-1:2*nod1, 2*nod2-1:2*nod2) = Mg(2*nod1-1:2*nod1, 2*nod2-1:2*nod2) + me(1:2, 3:4); 
	Mg(2*nod2-1:2*nod2, 2*nod1-1:2*nod1) = Mg(2*nod2-1:2*nod2, 2*nod1-1:2*nod1) + me(3:4, 1:2); 
	Mg(2*nod2-1:2*nod2, 2*nod2-1:2*nod2) = Mg(2*nod2-1:2*nod2, 2*nod2-1:2*nod2) + me(3:4, 3:4); 

end

% Construcao das matrizes para os elementos de viga horizontal.

for e=4:(Nel_vh-1)

	k0 = [12           6*data.L(e)     -12           6*data.L(e);
	      6*data.L(e)  4*(data.L(e))^2 -6*data.L(e)  2*(data.L(e))^2;
	      -12         -6*data.L(e)      12          -6*data.L(e);
	      6*data.L(e)  2*(data.L(e))^2 -6*data.L(e)  4*(data.L(e))^2];


	m0 = [156            22*data.L(e)      54           -13*data.L(e);
	      22*data.L(e)   4*(data.L(e))^2   13*data.L(e) -3*(data.L(e))^2;
	      54             13*data.L(e)      156          -22*data.L(e);
	      -13*data.L(e) -3*(data.L(e))^2  -22*data.L(e)  4*(data.L(e))^2];

	ke = data.E(e)*data.I(e)/(data.L(e)^3)*k0;
	me = (data.rho(e)*data.A(e)*data.L(e)/420)*m0;

	nod1 = con(e, 1);
	nod2 = con(e, 2);

	Kg(2*nod1-1:2*nod1, 2*nod1-1:2*nod1) = Kg(2*nod1-1:2*nod1, 2*nod1-1:2*nod1) + ke(1:2, 1:2); 
	Kg(2*nod1-1:2*nod1, 2*nod2-1:2*nod2) = Kg(2*nod1-1:2*nod1, 2*nod2-1:2*nod2) + ke(1:2, 3:4); 
	Kg(2*nod2-1:2*nod2, 2*nod1-1:2*nod1) = Kg(2*nod2-1:2*nod2, 2*nod1-1:2*nod1) + ke(3:4, 1:2); 
	Kg(2*nod2-1:2*nod2, 2*nod2-1:2*nod2) = Kg(2*nod2-1:2*nod2, 2*nod2-1:2*nod2) + ke(3:4, 3:4); 
	
	Mg(2*nod1-1:2*nod1, 2*nod1-1:2*nod1) = Mg(2*nod1-1:2*nod1, 2*nod1-1:2*nod1) + me(1:2, 1:2); 
	Mg(2*nod1-1:2*nod1, 2*nod2-1:2*nod2) = Mg(2*nod1-1:2*nod1, 2*nod2-1:2*nod2) + me(1:2, 3:4); 
	Mg(2*nod2-1:2*nod2, 2*nod1-1:2*nod1) = Mg(2*nod2-1:2*nod2, 2*nod1-1:2*nod1) + me(3:4, 1:2); 
	Mg(2*nod2-1:2*nod2, 2*nod2-1:2*nod2) = Mg(2*nod2-1:2*nod2, 2*nod2-1:2*nod2) + me(3:4, 3:4); 
end

% Construcao das matrizes para os elementos de portico (viga vertical).

% Preciso perguntar pro monitor pq nao sei se to fazendo certo ate aqui e
% de qualquer forma nao sei como seria a implementacao desse portico em EF.

% Passo 6 - Aplicar as condicoes de contorno para as fixacoes (na matriz de rigidez global).

Kgm = Kg; % Matriz de rigidez modificada.

% Basicamente zerar as linhas e colunas correspondentes aos
% nodes fixados. No caso aqui, sao ?.

%Kgm(:, 1) = 0; Kgm(1, :) = 0; Kgm(1, 1) = 1;
%Kgm(:, 2*1) = 0; Kgm(2*1, :) = 0; Kgm(2*1, 1) = 1;
%Kgm(:, 2*3-1) = 0; Kgm(2*3-1, :) = 0; Kgm(2*3-1, 2*3-1) = 1;
%Kgm(:, 2*3) = 0; Kgm(2*3, :) = 0; Kgm(2*3, 2*3) = 1;
