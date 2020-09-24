function [ output ] = bullseye( n )
mat = ones(2*n - 1);
a = size(mat,1);
b = size(mat,2);
c = n;
for i = 0:c
    for j = i + 1:a
        for k = i + 1:b
            mat(j,k) = c;
        end
    end
    a = a - 1;
    b = b - 1;
    c = c - 1;
end
output = mat;
end

    
    
    
            
    


