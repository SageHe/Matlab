clear all;
close all;
clc

opt_step = 1:101600;
res = .01/25.4; %mm/step
opt_translate = opt_step*res; %mm of translation

figure
hold on
plot(opt_step,opt_translate)
title ('Linear Translation Vs Steps','FontSize',18)
xlabel('Steps','fontsize',18)
ylabel('Linear Translation (in)','FontSize',18)
grid on
grid minor

dist = 0;
i = 1;
while i < 101600
    if(randi(100)==1) %
        dist(i) = dist(i-1);
    else
        dist(i) = dist(end)+(.01/25.4);
    end
    i = i + 1;
end
plot(dist)
legend('Optimal Translation','Translation W/skipped steps','FontSize',14)