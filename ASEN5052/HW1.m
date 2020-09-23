mu = 4e5;
R = 6378;
E = [1 100 0 0 10];
h = [1e5 1e5 7e4 8e4 8e4];
for i= 1:5
    if E(i) ~= 0
        rp(i) = (mu/(2*E(i)))*((sqrt(1 + ((2*E(i)*h(i)^2)/mu^2))) - 1);
    else
        rp(i) = .5*(h(i)^2/mu);
    end
    v(i) = sqrt((E(i) + (mu/R))*2);
end