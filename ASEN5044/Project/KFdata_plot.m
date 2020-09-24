clear;close all;clc
data = load('orbitdeterm_finalproj_KFdata.mat');
y = nan(36,1401);
for i = 2:size(data.ydata,2)
    vec = data.ydata{i};
    for j = 1:size(vec,2)
        y([vec(end,j)*3-2:vec(end,j)*3],i) = vec([1:3],j);
    end
end
figure
subplot(3,1,1)
for i=1:12
    scatter(data.tvec,y(3*i-2,:))
    hold on
end
subplot(3,1,2)
for i=1:12
    scatter(data.tvec,y(3*i-1,:))
    hold on
end
subplot(3,1,3)
for i=1:12
    scatter(data.tvec,y(3*i,:))
    hold on
end