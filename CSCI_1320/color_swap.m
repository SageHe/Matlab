function [out] = color_swap(current_img, r1, g1, b1, allowed, r2, g2, b2)
% [a,b,c] = size(current_img);
rd = (r1 - allowed:r1 + allowed);
gr = (g1 - allowed:g1 + allowed);
bl = (b1 - allowed:b1 + allowed);
out = current_img;
for i = 1:size(current_img,1)
    for j = 1:size(current_img,2)
        for k = 1:allowed*2
            if current_img(i,j,1) == rd(k)
                for L = 1:allowed*2
                    if current_img(i,j,2) == gr(L)
                        for m = 1:allowed*2
                            if current_img(i,j,3) == bl(m)
                                out(i,j,1) = r2;
                                out(i,j,2) = g2;
                                out(i,j,3) = b2;
                                %Using this series of nested for loops with
                                %embedded if statements, it doens't change
                                %the values specified by the user to the
                                %new values specified by the user unless
                                %all the the statements are met or made
                                %true in all of the nested loops and
                                %statements
                            end
                        end
                    end
                end
            end
        end
    end
end
end