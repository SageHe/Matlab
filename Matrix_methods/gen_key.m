clear all;close all;clc
temp = zeros(2);
keys = [];
for i = 0:25
    temp(1) = i;
    for j = 0:25 
        temp(3) = j;
        for k = 0:25
            temp(2) = k;
            for l = 0:25
                temp(4) = l;
                if ((gcd(round(det(temp)),13)) == 1) && ((gcd(round(det(temp)),2)) == 1) 
                    keys = cat(3,keys,temp);
                end
            end
        end
    end
end