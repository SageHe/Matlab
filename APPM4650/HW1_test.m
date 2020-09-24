%%
x = 1;
stuff = x*(x - 1)*(x - 2);

while stuff < 5.9994e15
    x = x+1;
    stuff = x*(x - 1)*(x - 2);
end

x = x - 1;
x
%% other question 
k = 1;
stuff = 1/((2*k + 1)*factorial(k));
stuff = stuff*(2/sqrt(pi));
while ((stuff > 1e-7) || (stuff < 0))
    k = k+1;
    stuff = (-1)^k/((2*k + 1)*factorial(k));
    stuff = stuff*(2/sqrt(pi));
end
clear stuff
stuff = 0;
for i = 0:10
    stuff = stuff + ((-1)^i*(1)^(2*i + 1))/((2*i + 1)*factorial(i));
end
stuff = stuff*(2/sqrt(pi));
stuff = 0;
denom = 1;
for i = 0:10
    denom = denom*(2*i + 1);
    stuff = stuff + (2^i)/denom;
end
stuff = stuff*(2/sqrt(pi))*exp(-1)
    
    
    