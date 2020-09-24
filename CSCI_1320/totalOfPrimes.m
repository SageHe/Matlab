function out = totalOfPrimes(m,n)
vec = m:n;
k = 1;
for i = m:1:n
    for j = 2:(i - 1)
            if mod(i,j) == 0 
                vec(k) = 0;
            end
    end
    k = k + 1;
end
out = sum(vec);
a = find(vec == 1);
vec(a) = [];
out = sum(vec);
end
    
            