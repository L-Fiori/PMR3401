clear; clc; close all;

b1 = 1.8; h1 = 0.9;
area1 = b1*h1;
b2 = 0.9; h2 = 0.9;
area2 = b2*h2;
D3 = 0.050;

dx = 1;

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
	Nod_p = size(0:dx:22, 2) + 2;
    
	coords = [0 5;
              2 5;
              4 5;
              5 5;
              6 5;
              (8:dx:64)' 5*ones(size(8:dx:64, 2), 1);
              65  5;
              36*ones(size(0:dx:22, 2), 1) (0:dx:22)';
              36 23];
          
	con = [ 4   48;
            11  48;
            35  48;
            (1:34)' (2:35)';
            36  37;
            37  38;
            38  20;
            20  39;
            (39:47)' (40:48)'];

elseif (dx == 3)
	Nod_vh = size(8:dx:63,2) + 5;
	Nod_p = size(0:dx:21, 2) + 2;
    
	coords = [ 0 5;
               3 5;
               5 5;
               6 5;
               (9:dx:63)' 5*ones(size(8:dx:63,2), 1);
               65 5;
               36*ones(size(0:dx:21, 2), 1) (0:dx:21)';
               36 23];              
          
	con = [ 4   33;
            8   33;
            24  33;
            (1:23)' (2:24)';
            25 26;
            26 14;
            14 27;
            (27:32)'  (28:33)'];

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

% Construcao das matrizes para os elementos de trelica (cabos)

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

% Construcao das matrizes para os elementos de portico (plataforma e torre).

for e=4:(4 + Nel_vh + Nel_p - 1)
    
    % Passo 3 - Construir a matriz de cada elemento.
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

% Passo 6 - Aplicar as condicoes de contorno para as fixacoes (na matriz de rigidez global).

Kgm = Kg; % Matriz de rigidez modificada.
Mgm = Mg; % Matriz de massa modificada.

% Basicamente zerar as linhas e colunas correspondentes aos
% nodes fixados. No caso aqui, sao todos os graus de liberdade do
% primeiro node da torre e translacoes do ultimo node da plataforma.

list = (1:Ngdl)';
i_cc = [(Nod_vh)*3 - 2, (Nod_vh)*3 - 1, (Nod_vh + 1)*3 - 2, (Nod_vh + 1)*3 - 1, (Nod_vh + 1)*3];

id_free = list(ismember(list, i_cc) == 0);
Kgm = Kg(id_free, id_free);
Mgm = Mg(id_free, id_free);

% Passo 7 - Encontrar os autovalores de [K -(omega^2)*M]

[A, V] = eigs(Mgm\Kgm, 6, 'smallestabs');
omega = sqrt(diag(V))/(2*pi)

% Passo 8 - Plotar os modos de vibracao

if (dx == 1)

    mod = 1
    U = A(:, mod);

    U = U';

    for i=1:5
        U = [U(1:i_cc(i)-1) 0 U(i_cc(i):end)];    
    end

    U = U';

    freq = omega(mod);
    T = 1/freq;
    dt = T/15;
    TF = 5*T;
    t = 0:dt:TF;
    scale = 10;
    
    coorExag = coords + scale*[ U(1:3:end) U(2:3:end) ]*sin(2*pi*freq*t(5));
    
    figure;
    plot_struct(coords, coorExag, con, '-r');
    
end

% =================== Analise Transiente ===================

a0 = 0.1217;
a1 = 0.0087;

% Retirando o terceiro grau de liberdade das matrizes pois nao serao utilizados.
Kg2 = Kgm;
Mg2 = Mgm;

% Matriz de amortecimento
Cgm = a0*Mg2 + a1*Kg2;

va = 1;

TF = 55/va;

dt_vec = [0.1; 0.01; 0.005];
dt1 = 0.1; dt2 = 0.01; dt3 = 0.005;
t1 = 0:dt1:TF;
t2 = 0:dt2:TF;
t3 = 0:dt3:TF;

V = zeros(26, 1);
Fnod = zeros(size(Mg2, 1), 1);

beta = 0.25;
gamma = 0.5;

% Prealocacao de vetor para deslocamento do node A para diferentes dts.
UA_dt1 = zeros(length(t1), 1);
UA_dt3 = zeros(length(t3), 1);

% Prealocacao de vetores para deslocamento dos nodes A, B, C e F para dt
% igual a 0,01.

UA = zeros(length(t2), 1);
UB = zeros(length(t2), 1);
UC = zeros(length(t2), 1);
UF = zeros(length(t2), 1);

% Loop para aplicacao dos carregamentos e implementacao do metodo Newmark
% Beta para aceleracao media constante.

for k=1:3
    if k == 1
        t = t1;
    elseif k == 2
        t = t2;
    elseif k == 3
        t = t3;
    end
    
    dt = dt_vec(k, 1)
    
    Meq = Mg2 + dt*gamma*Cgm + dt*dt*beta*Kg2;
    
    U1 = zeros(size(Mg2, 1), 1);
    V1 = zeros(size(Mg2, 1), 1);
    A1 = zeros(size(Mg2, 1), 1);
    
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

        Fnod(2:3:29) = -q1;
        Fnod(32:3:107) = V;
        Fnod(110:3:137) = -q2;

        F2 = Fnod;
        
        Feq = F2 - Cgm*( V1 + dt*(1-gamma)*A1 ) - Kg2*( U1 + dt*V1 + dt*dt*0.5*(1 - 2*beta)*A1 );
        
        A2 = Meq\Feq;
            
        U2 = U1 + dt*V1 + dt*dt*0.5*( (1-2*beta)*A1 + 2*beta*A2 );
        V2 = V1 + dt*( (1-gamma)*A1 + gamma*A2 );

        if k == 1
            UA_dt1(i, 1) = U2(2, 1);
        elseif k == 2
            UA(i, 1) = U2(2, 1);
            UB(i, 1) = U2(17, 1);
            UC(i, 1) = U2(29, 1);
            UF(i, 1) = U2(260, 1);
        elseif k == 3
            UA_dt3(i, 1) = U2(2, 1);
        end
            
        U1 = U2; V1 = V2; A1 = A2;
    end
    
end

figure;
plot(t2, UA)
hold on
plot(t2, UB)
hold on
plot(t2, UC)
hold on
plot(t2, UF)
legend('Ponto A', 'Ponto B', 'Ponto C', 'Ponto F')
title('Análise transiente para os pontos A, B, C e F, com passo dt = 0,01s')

figure;
plot(t1, UA_dt1)
hold on
plot(t2, UA)
hold on
plot(t3, UA_dt3)
legend('dt = 0.1s', 'dt = 0.01s', 'dt = 0.005s')
title('Análise transiente para o ponto A com diferentes passos')

% =================== Analise Harmonica ===================

df = 0.001;
freq = 0:df:6;
omega = 2*pi*freq;

F0 = zeros(size(Mg2, 1), 1);
F0(32:3:107) = -784;

YA = zeros(length(freq), 1);
YB = zeros(length(freq), 1);
YC = zeros(length(freq), 1);
YF = zeros(length(freq), 1);

for i=1:length(freq)
    Keq = (-omega(i)^2)*Mg2 + Kg2;
    
    Y = Keq\F0;
    
    YA(i, 1) = abs(Y(2, 1));
    YB(i, 1) = abs(Y(17, 1));
    YC(i, 1) = abs(Y(29, 1));
    YF(i, 1) = abs(Y(260, 1));
    
end

figure;
plot(freq, YA)
hold on
plot(freq, YB)
hold on
plot(freq, YC)
hold on
plot(freq, YF)
legend('Ponto A', 'Ponto B', 'Ponto C', 'Ponto F')
title('Diagrama de resposta em frequência para os pontos A, B, C e F')
