function [out] = truth(input)
for i = 1:length(input)
    if input(i,1) == 0
        if input(i,3) == 1
            s(i) = 1;
        else 
            s(i) = 0;
        end
    end
    if input(i,1) == 1
        if((input(i,2))||(input(i,3))) == 1
            s(i) = 1;
        elseif ((input(i,2))&&(input(i,3))) == 1
            s(i) = 1
        elseif ((input(i,2))&&(input(i,3))) == 0
            s(i) = 0;
        end
    end
end
out = s';
end
        