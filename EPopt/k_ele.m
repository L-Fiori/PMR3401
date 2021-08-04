function Ke = k_ele(mesh, element)

nodes = mesh.con(element, :);

x1 = mesh.coor(nodes(1), 1);
y1 = mesh.coor(nodes(1), 2);

x2 = mesh.coor(nodes(2), 1);
y2 = mesh.coor(nodes(2), 2);

x3 = mesh.coor(nodes(3), 1);
y3 = mesh.coor(nodes(3), 2);

% Calculo das constantes bi, ci e area Ae do elemento
b1 = y2 - y3;
b2 = y3 - y1;
b3 = y1 - y2;
c1 = x3 - x2;
c2 = x1 - x3;
c3 = x2 - x1;
Ae = (b1*c2 - b2*c1)/2;

% Calculo da matriz local do elemento
Ke = [(b1*b1 + c1*c1)  (b1*b2 + c1*c2)  (b1*b3 + c1*c3);
      (b1*b2 + c1*c2)  (b2*b2 + c2*c2)  (b2*b3 + c2*c3);
      (b1*b3 + c1*c3)  (b2*b3 + c2*c3)  (b3*b3 + c3*c3)];

Ke = (1/(4*Ae))*Ke;

end