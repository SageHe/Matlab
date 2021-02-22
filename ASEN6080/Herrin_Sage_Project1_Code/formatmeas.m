function [data] = formatmeas(meas)
data = NaN(meas(end,1)/20+1,size(meas,2));
data(:,1) = [0:20:meas(end,1)]';
data(1,:) = meas(1,:);
for i = 2:size(meas,1)
    data((meas(i,1)/20)+1,:) = meas(i,:);
end
end