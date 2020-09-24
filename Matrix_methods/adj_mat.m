%Takes any NXN matrix and determines its respective adjugate 
function [adj] = adj_mat(A)
    adj = zeros(size(A,1));
    %Preallocate adjugate matrix of zeroes
    for i = 1:size(A,1)
        for j = 1:size(A,1)
            temp = A;
            temp(j,:) = [];
            temp(:,i) = [];
            %'deletes' necessary row and column of original matrix
            B = det(temp);
            %Takes determinant of matrix 'temp'
            B = B*(-1)^(i + j);
            %Determines correct sign of adjugate entry based on 
            %entry position
            adj(i,j) = B;
            %Places B in i,j position in adjugate 
        end
    end
end
        
    
    