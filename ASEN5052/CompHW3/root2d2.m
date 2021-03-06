function F = root2d2(in)
mu = .5*(1 - sqrt(23/27)) + 1e-2;
x = in(1);
y = in(2);

F(1) = x + ((mu - 1).*(x+mu))./(((x+mu).^2 + y^2).^(3/2)) - (mu.*(x - 1 + mu))./(((x - 1 + mu).^2 + y^2).^(3/2));
F(2) = y + ((mu-1).*y)./(((x + mu).^2 + y^2).^(3/2)) - (mu.*y)./(((x - 1 + mu).^2 + y^2).^(3/2));
end