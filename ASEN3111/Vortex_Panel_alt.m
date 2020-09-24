function [c_l,x1,x2,x3,alpha,Cp,x] = Vortex_Panel_alt(x1,x2,x3,c,N,alpha)
[x,y,n] = NACA_Airfoils_alt(x1,x2,x3,c,N);
% y = [0 -.005 -.017 -.033 -.042 -.033 0 .045 .076 .072 .044 .013 0];
alpha = (alpha*pi)/180;
for i = 1:n
    xcp(i) = .5*(x(i) + x(i + 1));
    ycp(i) = .5*(y(i) + y(i + 1));
    panel_length(i) = sqrt((x(i + 1) - x(i))^2 + (y(i + 1) - y(i))^2);
    theta(i) = atan2((y(i + 1) - y(i)),(x(i + 1) - x(i)));
    rhs(i) = sin(theta(i) - alpha); %Rename this so not same as book
end
for i = 1:n
    for j = 1:n
        if i == j
            Cn_1(i,j) = -1;
            Cn_2(i,j) = 1;
            Ct_1(i,j) = pi/2;
            Ct_2(i,j) = pi/2;
        else
            A = -(xcp(i) - x(j))*cos(theta(j)) - (ycp(i) - y(j))*sin(theta(j));
            B = (xcp(i) - x(j))^2 + (ycp(i) - y(j))^2;
            C = sin(theta(i) - theta(j));
            D = cos(theta(i) - theta(j));
            E = (xcp(i) - x(j))*sin(theta(j)) - (ycp(i) - y(j))*cos(theta(j));
            F = log( 1 + panel_length(j)*(panel_length(j) + 2.*A)/B);
            G = atan2((E*panel_length(j)),B + A*panel_length(j));
            P = (xcp(i) - x(j))*sin(theta(i) - 2.*theta(j))...
                + (ycp(i) - y(j))*cos(theta(i) - 2.*theta(j));
            Q = (xcp(i) - x(j))*cos(theta(i) - 2.*theta(j))...
                - (ycp(i) - y(j))*sin(theta(i) - 2.*theta(j));
            Cn_2(i,j) = D + .5*Q*F/panel_length(j) - (A*C+D*E)*G/panel_length(j);
            Cn_1(i,j) = (.5*D*F) + C*G - Cn_2(i,j);
            Ct_2(i,j) = C + .5*P*F/panel_length(j) + (A*D-C*E)*G/panel_length(j);
            Ct_1(i,j) = .5*C*F - D*G - Ct_2(i,j);
        end
    end
end
%Influence coefficients
for i = 1:n
    AN(i,1) = Cn_1(i,1);
    AN(i,(n + 1)) = Cn_2(i,n);
    AT(i,1) = Ct_1(i,1);
    AT(i,(n + 1)) = Ct_2(i,n);
    for j = 2:n
        AN(i,j) = Cn_1(i,j) + Cn_2(i, (j - 1));
        AT(i,j) = Ct_1(i,j) + Ct_2(i,(j - 1));
    end
end
AN((n + 1),1) = 1;
AN((n + 1),(n + 1)) = 1;
for j = 2:n
    AN((n + 1),j) = 0;
end
rhs(n + 1) = 0;
rhs = rhs';
gamma = AN\rhs;
for i = 1:n
    V(i) = cos(theta(i) - alpha);
    for j = 1:(n + 1)
        V(i) = (V(i) + AT(i,j)*gamma(j));
        Cp(i) = 1 - V(i).^2;
    end
end
V = V';
Cp = Cp';
for i = 1:n
    ds(i) = sqrt((x(i + 1) - x(i))^2 + (y(i + 1) - y(i))^2);
end
ds = ds';
Gamma = sum(V.*ds);
c_l = 2*Gamma;
alpha = alpha*(180/pi);
% figure
% plot(x(1:end - 1),Cp)
% axis ij
% title(['Cp vs X/C NACA ',num2str(x1) num2str(x2) num2str(x3), ' at \alpha = ', num2str(alpha),char(176)])
% xlabel('Normalized Chord Length[m]')
% ylabel('Cp')
end