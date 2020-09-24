function [H] = calcH(X,i,j)
    Re = 6378;
    We = (2*pi)/86400;

    theta0 = (j - 1)*(pi/6);
    Xs = Re*cos(We*(i*10) + theta0);
    Xdots = -Re*sin(We*(i*10)+theta0)*We;
    Ys = Re*sin(We*(i*10) + theta0);
    Ydots = Re*cos(We*(i*10) + theta0)*We; 
    rho = sqrt((X(1)-Xs)^2 + (X(3)-Ys)^2);
    rhodot = ((X(1) - Xs)*(X(2)-Xdots) + (X(3)-Ys)*(X(4)-Ydots))/rho;
    num = ((X(1) - Xs)*(X(2)-Xdots) + (X(3)-Ys)*(X(4)-Ydots));
    phi = atan2((X(3)-Ys),(X(1)-Xs));
%     phi_check = atan2((X_gt(3,i)-Ys),(X_gt(1,i)-Xs));
    theta_i = atan2(Ys,Xs);
    H = [(X(1)-Xs)/rho 0 (X(3)-Ys)/rho 0;
            (rho*(X(2)-Xdots)-(num)*((X(1)-Xs)/rho))/(rho^2) (X(1)-Xs)/rho  (rho*(X(4)-Ydots)-(num)*((X(3)-Ys)/rho))/(rho^2) (X(3)-Ys)/rho;
            (Ys-X(3))/(rho^2) 0 (X(1)-Xs)/(rho^2) 0];
end