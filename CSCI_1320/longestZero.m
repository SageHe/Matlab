function out = longestZero(S)
temp = [];
%counter = 0;
delete = [];
%maxiumum = 0;
for i = 1:length(S)
    if S(i) == '0'
        temp = [temp i];
        %counter = counter + 1;
        %max = counter;
    elseif S(i) == '1'
        if length(temp) > length(delete)
            delete = temp;
            temp = [];
            %counter = 0;
        else temp = [];
            %counter = 0;
        end
    end
end
if length(temp) > length(delete)
    delete = temp;
end
S(delete) = [];
out = S;
end
        
