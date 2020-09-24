R1 = [1 0 0;0 cosd(-60) sind(-60);0 -sind(-60) cosd(-60)];
R2 = [cosd(-60) 0 -sind(-60); 0 1 0; sind(-60) 0 cosd(-60)];
R3 = [cosd(-60) sind(-60) 0; -sind(-60) cosd(-60) 0; 0 0 1];

a = R3*R2*R1;
b = R1*R2*R3;
c = R3*R1*R3;

[V, D] = eig(c);

V

D
