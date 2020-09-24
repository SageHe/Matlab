function [out,removed_rows] = remove_rows(M)
removed_rows = [];
for i = 1:size(M,1)
    for j = 1:size(M,2)
        if M(i,j) == inf
            M(i,:) = zeros(1,size(M,2));
        end
        if isnan(M(i,j)) == 1
            M(i,:) = zeros(1,size(M,2));
        end
%         if i == (size(M,1) + 1)
%             return;
%         end
    end
end
for i = 1:size(M,1)
    if mean(M(i,:)) == 0
        removed_rows = [removed_rows i];
    end
end
M(all(M==0,2),:)=[];
out = M;
end