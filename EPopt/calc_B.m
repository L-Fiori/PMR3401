function B = calc_B(x, y)

x1 = x(1); x2 = x(2); x3 = x(3);
y1 = y(1); y2 = y(2); y3 = y(3);

b1 = y2 - y3;
b2 = y3 - y1;
c1 = x3 - x2;
c2 = x1 - x3;
Ae = (b1*c2 - b2*c1)/2;

dN1dx = (1/(2*Ae))*(y2 - y3);
dN2dx = (1/(2*Ae))*(y3 - y1);
dN3dx = (1/(2*Ae))*(y1 - y2);

dN1dy = (1/(2*Ae))*(x3 - x2);
dN2dy = (1/(2*Ae))*(x1 - x3);
dN3dy = (1/(2*Ae))*(x2 - x1);

B = [dN1dx dN2dx dN3dx;
     dN1dy dN2dy dN3dy];

end

