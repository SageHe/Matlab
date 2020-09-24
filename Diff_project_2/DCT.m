function out = DCT(mat)
DCT = [];
for i = 1:size(mat,1)
    for j = 1:size(mat,2)
        DCT(i,j) = sqrt((2/length(mat)))*cos((pi*(i - (1/2))*(j - (1/2)))/length(mat));
    end
end
out = DCT;
end
