clear; clc; close all;

b1 = 1.8; h1 = 0.9;
area1 = b1*h1;
b2 = 0.9; h2 = 0.9;
area2 = b2*h2;
D3 = 0.050;

dx = 2;

% Passo 1 - Definir a geometria (coordenadas e conectividades).

if (dx == 1)
	Nod_vh = size(0:dx:65, 2);
	Nod_p = size(0:dx:23, 2);
    
	coords = [(0:dx:65)' 5*ones(Nod_vh, 1);
              36*ones(Nod_p, 1) (0:dx:23)'];
          
    coords(72, :) = [];
    
    con = [ 6 89;
        19 89;
        66 89;
        (1:(Nod_vh - 1))' (2:(Nod_vh))';
        (67:70)' (68:71)';
        71 37;
        37 72;
        (72:88)' (73:89)'];
    
elseif (dx == 2)
	Nod_vh = size(8:dx:64, 2) + 6;
	Nod_p = size(0:dx:22, 2) + 1;
    
	coords = [0 5
              2 5
              4 5
              5 5
              6 5
              (8:dx:64)' 5*ones(size(8:dx:64, 2), 1);
              65  5;
              36*ones(size(0:dx:22, 2), 1) (0:dx:22)';
              36 23];
          
	con = [ 4   48;
            11  48;
            35  48;
            (1:34)' (2:35)';
            (36:47)' (37:48)'];

elseif (dx == 3)
	Nod_vh = size(0:dx:65, 2) + 1;
	Nod_p = size(0:dx:23, 2) + 1;
    
	coords = [(0:dx:65)' 5*ones(Nod_vh-1, 1);
              65  5;
              36*ones(Nod_p-1, 1) (0:dx:23)';
              36 23
              5 5];
          
	con = [ 33  32;
            7  32;
            23  32;
            (1:22)' (2:23)';
            (24:31)' (25:32)'];

end
        
%plot_struct(coords, con, '-b');


% Passo 2 - Atribuir propriedades aos elementos.

Nod = size(coords, 1);
Nel = size(con, 1);
Nel_vh = Nod_vh - 1;
Nel_p = Nod_p - 1;
Ngdl = 3*Nod;

data.E = 210e9*ones(Nel, 1); % Vetor de modulos de Young dos elementos.
data.rho = 7600*ones(Nel, 1); % Vetor de massa especifica dos elementos.
data.A = ones(Nel, 1); % Vetor de areas dos elementos.
data.I = ones(Nel, 1); % Vetor de momentos de inercia dos elementos.

data.A(1:3) = pi*(D3/2)^2;
data.A(4:(4 + Nel_vh - 1)) = area2;
data.A((4 + Nel_vh):(4 + Nel_vh + Nel_p - 1)) = area1;

data.I(4:(4 + Nel_vh - 1)) = b2*((h2)^3)/12;
data.I((4 + Nel_vh):(4 + Nel_vh + Nel_p - 1)) = h1*((b1)^3)/12;

data.L = zeros(Nel, 1); % Vetor de comprimento dos elementos.
data.Q = zeros(Nel, 1); % Vetor de rotacoes dos elementos.

for e=1:Nel

	x1 = coords(con(e, 1), 1);
	x2 = coords(con(e, 2), 1);

	y1 = coords(con(e, 1), 2);
	y2 = coords(con(e, 2), 2);

	data.L(e) = sqrt( (x2 - x1)^2 + (y2 - y1)^2 );
	data.Q(e) = atan2((y2 - y1), (x2 - x1))*(180/pi);

end

% Matrizes globais de rigidez e massa.

Kg = zeros(Ngdl, Ngdl);
Mg = zeros(Ngdl, Ngdl);

% Construcao das matrizes para os elementos de trelica

k0_t = [ 1 0 0 -1 0 0;
         0 0 0  0 0 0;
	     0 0 0  0 0 0;
	    -1 0 0  1 0 0;
	     0 0 0  0 0 0;
	     0 0 0  0 0 0];

m0_t = [2 0 0 1 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0;
        1 0 0 2 0 0;
        0 0 0 0 0 0;
        0 0 0 0 0 0];
    

for e=1:3

	% Passo 3 - Construir a matriz de cada elemento.
	ke = (data.E(e)*data.A(e)/data.L(e))*k0_t;
	me = (data.rho(e)*data.A(e)*data.L(e)/6)*m0_t;

	% Passo 4 - Rotacionar cada elemento para o sistema global.
	c = cosd(data.Q(e)); s = sind(data.Q(e));
	T = [ c s 0  0 0 0;
	     -s c 0  0 0 0;
	      0 0 1  0 0 0;
	      0 0 0  c s 0;
	      0 0 0 -s c 0;
	      0 0 0  0 0 1];
	ke = T'*ke*T;
	me = T'*me*T;

	% Passo 5 - Alocar a informacao de cada matriz de elemento na matriz global.
	nod1 = con(e, 1);
	nod2 = con(e, 2);

	Kg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) = Kg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) + ke(1:3, 1:3); 
	Kg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) = Kg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) + ke(1:3, 4:6); 
	Kg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) = Kg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) + ke(4:6, 1:3); 
	Kg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) = Kg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) + ke(4:6, 4:6); 
	
	Mg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) = Mg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) + me(1:3, 1:3); 
	Mg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) = Mg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) + me(1:3, 4:6); 
	Mg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) = Mg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) + me(4:6, 1:3); 
	Mg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) = Mg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) + me(4:6, 4:6); 

end

% Construcao das matrizes para os elementos de viga horizontal.
%{
for e=4:(4 + Nel_vh - 1)

	k0_v = [0  0           0               0  0            0;
		    0  12          6*data.L(e)     0 -12           6*data.L(e);
	        0  6*data.L(e) 4*(data.L(e))^2 0 -6*data.L(e)  2*(data.L(e))^2;
	        0  0           0               0  0            0;       
		    0 -12         -6*data.L(e)     0  12          -6*data.L(e);
	        0  6*data.L(e) 2*(data.L(e))^2 0 -6*data.L(e)  4*(data.L(e))^2];

	m0_v = [ 0  0              0                0 0              0;
		     0  156            22*data.L(e)     0 54            -13*data.L(e);
	         0  22*data.L(e)   4*(data.L(e))^2  0 13*data.L(e)  -3*(data.L(e))^2;
	         0  0              0                0 0              0;       
		     0  54             13*data.L(e)     0 156           -22*data.L(e);
	         0 -13*data.L(e)  -3*(data.L(e))^2  0 -22*data.L(e)  4*(data.L(e))^2];

	ke = data.E(e)*data.I(e)/(data.L(e)^3)*k0_v;
	me = (data.rho(e)*data.A(e)*data.L(e)/420)*m0_v;

	nod1 = con(e, 1);
	nod2 = con(e, 2);

	Kg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) = Kg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) + ke(1:3, 1:3); 
	Kg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) = Kg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) + ke(1:3, 4:6); 
	Kg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) = Kg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) + ke(4:6, 1:3); 
	Kg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) = Kg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) + ke(4:6, 4:6); 
	
	Mg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) = Mg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) + me(1:3, 1:3); 
	Mg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) = Mg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) + me(1:3, 4:6); 
	Mg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) = Mg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) + me(4:6, 1:3); 
	Mg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) = Mg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) + me(4:6, 4:6); 

end
%}
% Construcao das matrizes para os elementos de portico (viga vertical).

        
%for e=(4 + Nel_vh):(4 + Nel_vh + Nel_p - 1)

for e=4:(4 + Nel_vh + Nel_p - 1)
	ke_t = (data.E(e)*data.A(e)/data.L(e))*k0_t;
	me_t = (data.rho(e)*data.A(e)*data.L(e)/6)*m0_t;

	k0_v = [0  0           0               0  0            0;
		    0  12          6*data.L(e)     0 -12           6*data.L(e);
	        0  6*data.L(e) 4*(data.L(e))^2 0 -6*data.L(e)  2*(data.L(e))^2;
	        0  0           0               0  0            0;       
		    0 -12         -6*data.L(e)     0  12          -6*data.L(e);
	        0  6*data.L(e) 2*(data.L(e))^2 0 -6*data.L(e)  4*(data.L(e))^2];

	m0_v = [ 0  0              0                0 0              0;
		     0  156            22*data.L(e)     0 54            -13*data.L(e);
	         0  22*data.L(e)   4*(data.L(e))^2  0 13*data.L(e)  -3*(data.L(e))^2;
	         0  0              0                0 0              0;       
		     0  54             13*data.L(e)     0 156           -22*data.L(e);
	         0 -13*data.L(e)  -3*(data.L(e))^2  0 -22*data.L(e)  4*(data.L(e))^2];

	ke_v = (data.E(e)*data.I(e)/(data.L(e)^3))*k0_v;
	me_v = (data.rho(e)*data.A(e)*data.L(e)/420)*m0_v;

	ke = ke_t + ke_v;
	me = me_t + me_v;

	c = cosd(data.Q(e)); s = sind(data.Q(e));
	T = [ c s 0  0 0 0;
	     -s c 0  0 0 0;
	      0 0 1  0 0 0;
	      0 0 0  c s 0;
	      0 0 0 -s c 0;
	      0 0 0  0 0 1];

	ke = T'*ke*T;
	me = T'*me*T;

	nod1 = con(e, 1);
	nod2 = con(e, 2);
    
    Kg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) = Kg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) + ke(1:3, 1:3); 
	Kg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) = Kg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) + ke(1:3, 4:6); 
	Kg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) = Kg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) + ke(4:6, 1:3); 
	Kg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) = Kg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) + ke(4:6, 4:6); 
	
	Mg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) = Mg(3*nod1-2:3*nod1, 3*nod1-2:3*nod1) + me(1:3, 1:3); 
	Mg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) = Mg(3*nod1-2:3*nod1, 3*nod2-2:3*nod2) + me(1:3, 4:6); 
	Mg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) = Mg(3*nod2-2:3*nod2, 3*nod1-2:3*nod1) + me(4:6, 1:3); 
	Mg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) = Mg(3*nod2-2:3*nod2, 3*nod2-2:3*nod2) + me(4:6, 4:6); 
    
    %{
	Kg((2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln)), (2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln))) = Kg((2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln)), (2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln))) + ke(1:3, 1:3); 
	Kg((2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln)), (2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln))) = Kg((2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln)), (2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln))) + ke(1:3, 4:6); 
	Kg((2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln)), (2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln))) = Kg((2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln)), (2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln))) + ke(4:6, 1:3); 
	Kg((2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln)), (2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln))) = Kg((2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln)), (2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln))) + ke(4:6, 4:6); 
	
	Mg((2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln)), (2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln))) = Mg((2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln)), (2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln))) + me(1:3, 1:3); 
	Mg((2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln)), (2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln))) = Mg((2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln)), (2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln))) + me(1:3, 4:6); 
	Mg((2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln)), (2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln))) = Mg((2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln)), (2*v_ln + 3*(nod1-v_ln) - 2):(2*v_ln + 3*(nod1-v_ln))) + me(4:6, 1:3); 
	Mg((2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln)), (2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln))) = Mg((2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln)), (2*v_ln + 3*(nod2-v_ln) - 2):(2*v_ln + 3*(nod2-v_ln))) + me(4:6, 4:6); 
    %}
    
end

% Passo 6 - Aplicar as condicoes de contorno para as fixacoes (na matriz de rigidez global).

Kgm = Kg; % Matriz de rigidez modificada.
Mgm = Mg; % Matriz de massa modificada.

% Basicamente zerar as linhas e colunas correspondentes aos
% nodes fixados. No caso aqui, sao todos os graus de liberdade do
% node 5 + Nod_vh e w2 do node 5 + Nod_vh - 1 (ultimo node da viga
% horizontal).

%{
Kgm(:, (Nod_vh + 1)*3 - 2) = 0; Kgm((Nod_vh + 1)*3 - 2, :) = 0; Kgm((Nod_vh + 1)*3 - 2, (Nod_vh + 1)*3 - 2) = 1;
Kgm(:, (Nod_vh + 1)*3 - 1) = 0; Kgm((Nod_vh + 1)*3 - 1, :) = 0; Kgm((Nod_vh +1)*3 - 1, (Nod_vh + 1)*3 - 1) = 1;
Kgm(:, (Nod_vh + 1)*3) = 0; Kgm((Nod_vh + 1)*3, :) = 0; Kgm((Nod_vh + 1)*3, (Nod_vh + 1)*3) = 1;

Kgm(:, (Nod_vh)*3 - 2) = 0; Kgm((Nod_vh)*3 - 2, :) = 0; Kgm((Nod_vh)*3 - 2, (Nod_vh)*3 - 2) = 1;
Kgm(:, (Nod_vh)*3 - 1) = 0; Kgm((Nod_vh)*3 - 1, :) = 0; Kgm((Nod_vh)*3 - 1, (Nod_vh)*3 - 1) = 1;

Mgm(:, (Nod_vh + 1)*3 - 2) = 0; Mgm((Nod_vh + 1)*3 - 2, :) = 0; Mgm((Nod_vh + 1)*3 - 2, (Nod_vh + 1)*3 - 2) = 1;
Mgm(:, (Nod_vh + 1)*3 - 1) = 0; Mgm((Nod_vh + 1)*3 - 1, :) = 0; Mgm((Nod_vh +1)*3 - 1, (Nod_vh + 1)*3 - 1) = 1;
Mgm(:, (Nod_vh + 1)*3) = 0; Mgm((Nod_vh + 1)*3, :) = 0; Mgm((Nod_vh + 1)*3, (Nod_vh + 1)*3) = 1;

Mgm(:, (Nod_vh)*3 - 2) = 0; Mgm((Nod_vh)*3 - 2, :) = 0; Mgm((Nod_vh)*3 - 2, (Nod_vh)*3 - 2) = 1;
Mgm(:, (Nod_vh)*3 - 1) = 0; Mgm((Nod_vh)*3 - 1, :) = 0; Mgm((Nod_vh)*3 - 1, (Nod_vh)*3 - 1) = 1;

% Passo 7 - Encontrar os autovalores de [K -(omega^2)*M]

Kgm(:, (Nod_vh + 1)*3 - 2) = []; Kgm((Nod_vh + 1)*3 - 2, :) = [];
Kgm(:, (Nod_vh + 1)*3 - 1) = []; Kgm((Nod_vh + 1)*3 - 1, :) = [];
Kgm(:, (Nod_vh + 1)*3) = []; Kgm((Nod_vh + 1)*3, :) = [];

Kgm(:, (Nod_vh)*3 - 2) = []; Kgm((Nod_vh)*3 - 2, :) = [];
Kgm(:, (Nod_vh)*3 - 1) = []; Kgm((Nod_vh)*3 - 1, :) = [];

Mgm(:, (Nod_vh + 1)*3 - 2) = []; Mgm((Nod_vh + 1)*3 - 2, :) = [];
Mgm(:, (Nod_vh + 1)*3 - 1) = []; Mgm((Nod_vh + 1)*3 - 1, :) = [];
Mgm(:, (Nod_vh + 1)*3) = []; Mgm((Nod_vh + 1)*3, :) = [];

Mgm(:, (Nod_vh)*3 - 2) = []; Mgm((Nod_vh)*3 - 2, :) = [];
Mgm(:, (Nod_vh)*3 - 1) = []; Mgm((Nod_vh)*3 - 1, :) = [];
%}

i_cc = [(Nod_vh)*3 - 2, (Nod_vh)*3 - 1, (Nod_vh + 1)*3 - 2, (Nod_vh + 1)*3 - 1, (Nod_vh + 1)*3];

indices = []; % indices das colunas nulas - a serem removidas

%{
for i=1:size(Kgm, 1)
    if Kgm(i, :) == zeros(1, size(Kgm, 2))
        indices(end+1) = i;
    end
end
%}

%xlswrite('Kgm_zeros.xlsx', Mgm)

%Kgm(~any(Kgm, 2), :) = []; %rows
%Kgm(:, ~any(Kgm, 1)) = []; %columns

%xlswrite('Kgm.xlsx', Kgm)

%xlswrite('Mgm_zeros.xlsx', Mgm)

%Mgm(~any(Mgm, 2), :) = []; %rows
%Mgm(:, ~any(Mgm, 1)) = []; %columns

%xlswrite('Mgm.xlsx', Mgm)

% Exclui linhas que estão fixadas

list = (1:Ngdl)';
i_cc = [(Nod_vh)*3 - 2, (Nod_vh)*3 - 1, (Nod_vh + 1)*3 - 2, (Nod_vh + 1)*3 - 1, (Nod_vh + 1)*3];

id_fix = [((Nod_vh)*3 - 2) ((Nod_vh)*3-1) ((Nod_vh+1)*3-2) ((Nod_vh+1)*3-1) ((Nod_vh+1)*3)]';
id_free = list(ismember(list,id_fix) == 0);
Kgm = Kg(id_free,id_free);
Mgm = Mg(id_free,id_free);
%[A, V] = eigs(Mgm\Kgm, 6, 'smallestabs');

%omega_1 = sqrt(diag(V))/(2*pi)

A = Mgm\Kgm;
[vec, val] = eigs(A, 6, 'smallestabs'); %smallestabs
val = sqrt(diag(val))/(2*pi);
val

if (dx == 2)

    mod = 3
    U = A(:, mod);

    U = U';
    for i=2:size(indices, 2)
        U = [U(1:indices(i)-1) 0 U(indices(i):end)];
    end

    for i=1:5
        U = [U(1:i_cc(i)-1) 0 U(i_cc(i):end)];    
    end

    U = U';

    freq = omega_1(mod);
    T = 1/freq;
    dt = T/15;
    TF = 5*T;
    t = 0:dt:TF;
    scale = 10;

    for i=1:length(t)
        coorExag = coords + scale*[ U(1:3:end) U(2:3:end) ]*sin(2*pi*freq*t(i));
        plot_struct(coorExag, con, '-r');
        axis([-1 66 -1 24])
        pause(0.1);
        clf;
    end
    
end


a0 = 0.1217;
a1 = 0.0087;

Kg2 = Kgm;
Mg2 = Mgm; % Retirando o terceiro grau de liberdade das matrizes pois nao serao utilizadas.

%Kg2(3:3:end, :) = []; Kg2(:, 3:3:end) = [];
%Mg2(3:3:end, :) = []; Mg2(:, 3:3:end) = [];
%Mg2(131:134, :) = []; Mg2(:, 131:134) = [];
%Kg2(131:134, :) = []; Kg2(:, 131:134) = [];

%{
Kg2(~any(Kg2, 2), :) = []; %rows
Kg2(:, ~any(Kg2, 1)) = []; %columns

Mg2(~any(Mg2, 2), :) = []; %rows
Mg2(:, ~any(Mg2, 1)) = []; %columns
%}

%{
for i=1:size(Kg2, 1)
    if Kg2(i, :) == zeros(1, size(Kg2, 1))
        Kg2(i, i) = 1;
        Mg2(i, i) = 1;
    end
end
%}

Cgm = a0*Mg2 + a1*Kg2;

va = 2;

TF = 55/va;

dt = 0.05;
t = 0:dt:TF;

U1 = zeros(size(Mg2, 1), 1);
V1 = zeros(size(Mg2, 1), 1);
A1 = zeros(size(Mg2, 1), 1);

V = zeros(26, 1);
Fnod = zeros(size(Mg2, 1), 1);

beta = 0.25;
gamma = 0.5;

Meq = Mg2 + dt*gamma*Cgm + dt*dt*beta*Kg2;

UA = zeros(length(t), 1);
UB = zeros(length(t), 1);
UC = zeros(length(t), 1);
UF = zeros(length(t), 1);

for i=1:length(t)
    
    if (t(i) <= 27/va)
        N1 = 20;
    elseif t(i) <= 47/va
        N1 = va*t(i) - 7;
    elseif (t(i) > 47/va)
        N1 = 40;
    end
    
    if (t(i) <= 20/va)
        N2 = 20 - va*t(i);
    elseif (t(i) > 20/va)
        N2 = 0;
    end
    
    for j=1:size(V, 1)
    
        if ( t(i) <= (36 - (9 + j))/va )
            V(j, 1) = 0;
        elseif ( t(i) <= (36 + 20 -(9 + j))/va )
            V(j, 1) = -80*9.8*(1 - cos(2*pi*va*t(i)))/2;
        elseif (t(i) > (36 + 20 -(9 + j))/va)
            V(j, 1) = 0;
        end
        
    end
  
    q1 = N1*80*9.8/9;
    q2 = N2*80*9.8/9;
    
    %{
    Fnod(2:2:20) = -q1;
    Fnod(22:2:72) = V;
    Fnod(74:2:92) = -q2;
    %}
    
    Fnod(2:3:29) = -q1;
    Fnod(32:3:107) = V;
    Fnod(110:3:137) = -q2;
    
    F2 = Fnod;
    
    %{
    if (i == 1)
        U1 = Kg2\Fnod
    end
    %}
    
    Feq = F2 - Cgm*( V1 + dt*(1-gamma)*A1 ) - Kg2*( U1 + dt*V1 + dt*dt*0.5*(1 - 2*beta)*A1 );
    A2 = Meq\Feq;
    U2 = U1 + dt*V1 + dt*dt*0.5*( (1-2*beta)*A1 + 2*beta*A2 );
    V2 = V1 + dt*( (1-gamma)*A1 + gamma*A2 );
    
    
    UA(i, 1) = U2(2, 1);
    UB(i, 1) = U2(17, 1);
    UC(i, 1) = U2(29, 1);
    UF(i, 1) = U2(260, 1);
    
    %{
    UA(i, 1) = U2(1, 1);
    UB(i, 1) = U2(7, 1);
    UC(i, 1) = U2(11, 1);
    UF(i, 1) = U2(111, 1);
    %}
   
    U1 = U2; V1 = V2; A1 = A2;
    
end

%{
plot(t, UA)
hold on
plot(t, UB)
hold on
plot(t, UC)
hold on
plot(t, UF)
%}

