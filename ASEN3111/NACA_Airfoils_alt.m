%Author: Sage Herrin
%Created: 10/12/18
%SID:106071909
%Function that plots NACA airfoils based on 4 digit input. x1 is thickness
%distribution, given as 100m => m = x1/100, x2 is location of max camber p,
%s.t x2 = 10p => p = x2/10. x3 represents max thickness as function
%of fraction of chord s.t x3 = 100t => t = x3/100
function [x,y,n] = NACA_Airfoils_alt(x1,x2,x3,c,N)
%Define NACA defined variables
m = x1/100;
p = x2/10;
t = x3/100;
%Thickness distribution from mean camber
theta = linspace(-2*pi,0,(N + 1));
x = (c/2)*(cos(theta) + 1);
x_cam = linspace(0,c,(N + 1));
yt = (t/.2)*c*(.2969*sqrt(x./c) - .1260*(x./c) - .3516*(x./c).^2 + .2843*(x./c).^3 - .1036*(x./c).^4);

yc = [];
for i = 1:(N + 1)
    if x(i) < p*c
        yc = [yc m*(x(i)/(p^2))*(2*p - (x(i)/c))];
    else
        yc = [yc m*((c - x(i))/((1 - p)^2))*(1 + (x(i)/c) - 2*p)];
    end
end

xi = atan(gradient(yc)./gradient(x));
xi(find(isnan(xi))) = 0;
a = find(x == 0);
xi_upper = xi(1:a);
xi_lower = xi(a:end);

xu = x(1:a) - (yt(1:a).*sin(xi_upper));
yu = yc(1:a) + (yt(1:a).*cos(xi_upper));
xl = x(a:end) + (yt(a:end).*sin(xi_lower));
yl = yc(a:end) - (yt(a:end).*cos(xi_lower));

% xu = flip(xu);
% yl = flip(yl);
xl(1) = [];
yl(1) = [];

% x = [xu xl];
y = [yu yl];
y = flip(y);
% y(find(isnan(y))) = 0;
% figure
% hold on
% % plot(x_cam,yc)
% plot(x,yc)
% plot(x,y)
% % plot(xu,yu)
% % plot(xl,yl)
% axis 'equal'
% title(['NACA ', num2str(x1) num2str(x2) num2str(x3)])
% xlabel('X[m]')
% ylabel('Y[m]')
% legend('Mean Camber Line')
% hold off
n = N;
end

