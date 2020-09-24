theta = acos((-2*sqrt(2)/3));
W1 = [0 0 1];
W2 = [cos(theta) 0 -sin(theta)];
W3 = [cos(theta)*cosd(120) cos(theta)*sind(120) -sin(theta)];
W4 = [cos(theta)*cosd(120) -cos(theta)*sind(120) -sin(theta)];

angle1 = acosd((dot(W1,W2))/(norm(W1)*norm(W2)));
angle2 = acosd((dot(W2,W3))/(norm(W2)*norm(W3)));
angle3 = acosd((dot(W3,W4))/(norm(W3)*norm(W4)));
angle4 = acosd((dot(W4,W1))/(norm(W4)*norm(W1)));


A = [W2' W3' W4'];
b = [1 0 0]';

omega1 = A\b;

clear b;

b = [0 1 0]';

omega2 = A\b;