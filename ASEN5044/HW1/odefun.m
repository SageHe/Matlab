function dydt = odefun(t,y)
%Define w_n, equal to sqrt(k/m)
w_n = 1;
%Preallocate function dydt size
dydt = zeros(2,1);
%Define dydt first and second indeces 
dydt(1) = y(2);
dydt(2) = -y(1);
