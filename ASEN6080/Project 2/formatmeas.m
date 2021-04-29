function [data] = formatmeas(meas)
data(1,:) = meas(1,:);
ind = 1;
cnt = 0;
while ind < size(meas,1)
    if ((meas(ind+1,1)-data(end,1)) > 20) && ((meas(ind+1,1)-data(end,1)) > 1200)
        data(end+1,1) = data(end,1) + 1200;
        data(end,2:7) = NaN;
    else
        data(end+1,:) = meas(ind+1,:);
        ind = ind + 1;
    end
end 
end