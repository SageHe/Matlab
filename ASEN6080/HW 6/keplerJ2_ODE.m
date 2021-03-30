function zd = keplerJ2_ODE(t,Z,n)
Z = reshape(Z,n,2*n+1);
Z = Z';
% mu = 3.986004415e5;
% J2 = 1082.63e-6;
for i = 1:(2*n+1)
    x = Z(i,1);
    y = Z(i,2);
    z = Z(i,3);
    xdot = Z(i,4);
    ydot = Z(i,5);
    zdot = Z(i,6);
    mu = Z(i,7)*1e4;
    J2 = Z(i,8)*1e-3;
    rdot_vec = [xdot;ydot;zdot];

    r = sqrt(x^2 + y^2 + z^2);
    R = 6378.0;

    xdd = (mu*((J2*R^2*x)/(x^2 + y^2 + z^2)^2 - (2*J2*R^2*x*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (mu*x*((J2*R^2*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^2) + 1))/(x^2 + y^2 + z^2)^(3/2);
    ydd = (mu*((J2*R^2*y)/(x^2 + y^2 + z^2)^2 - (2*J2*R^2*y*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (mu*y*((J2*R^2*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^2) + 1))/(x^2 + y^2 + z^2)^(3/2);
    zdd = - (mu*((2*J2*R^2*z)/(x^2 + y^2 + z^2)^2 + (2*J2*R^2*z*(x^2 + y^2 - 2*z^2))/(x^2 + y^2 + z^2)^3))/(x^2 + y^2 + z^2)^(1/2) - (mu*z*((J2*R^2*(x^2 + y^2 - 2*z^2))/(2*(x^2 + y^2 + z^2)^2) + 1))/(x^2 + y^2 + z^2)^(3/2);

    vdot_vec = [xdd;ydd;zdd];

    xdot_vec = [rdot_vec;vdot_vec];

    zd(i,:) = [xdot_vec' 0 0];
end
zd = zd';
zd = reshape(zd,2*n^2 + n,1); 
end