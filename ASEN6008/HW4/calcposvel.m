function [R,V] = calcposvel(a,e,i,omega,womega,nu,mu);
% mu = 3.986004415e5;
% mu = 1.32712440018e11;
ehat = [(cosd(womega)*cosd(omega))-(cosd(i)*sind(womega)*sind(omega));...
        (cosd(womega)*sind(omega))+(cosd(i)*sind(womega)*cosd(omega));...
        (sind(womega)*sind(i))];
ehat_perp = [-((sind(womega)*cosd(omega))+(cosd(i)*cosd(womega)*sind(omega)));...
                -((sind(womega)*sind(omega))-(cosd(i)*cosd(womega)*cosd(omega)));...
                (cosd(womega)*sind(i))];
r = (a*(1 - e^2))/(1 + e*cosd(nu)); 

R = r*(cosd(nu)*ehat + sind(nu)*ehat_perp);

V = sqrt(mu/(a*(1-e^2)))*(-sind(nu)*ehat + (e + cosd(nu))*ehat_perp);
end