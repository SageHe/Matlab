function stat_num = DST(Y,ydata_statnum)
temp = ~isnan(Y);
if sum(temp) == size(ydata_statnum)*3
    stat_num = ydata_statnum;
    return
end
ind = find(temp == 1);
for i = 1:size(ind)/3
    Y_statnum(i) = ind(i*3)/3;
end
stat_num = intersect(Y_statnum,ydata_statnum);
end