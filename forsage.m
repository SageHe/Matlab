t = linspace(-pi,pi, 350);
X = t .* sin( pi * .872*sin(t)./t);
Y = -abs(t) .* cos(pi * sin(t)./t);
plot(X,Y);
fill(X, Y, 'r');
axis square;
title('hi', 'FontSize', 28);