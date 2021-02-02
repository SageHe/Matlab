function [out] = parsemeas(meas,t)
mat = [t meas];
meas = [];
for i = 1:numel(t)
    if (~isnan(mat(i,2))||~isnan(mat(i,3))||~isnan(mat(i,4)))
        meas = [meas;mat(i,:)];
    end
end
out = meas;
counter = 0;
for i = 1:size(meas,1)
    for j = 2:4
        if isnan(meas(i,j)) ~= 1
            counter = counter + 1;
        end
    end
end

end