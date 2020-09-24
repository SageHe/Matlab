x = 1:10;
y1 = x;
y2 = 2*x;
y3 = 3*x;
y4 = 4*x;
y5 = 5*x;
y6 = 6*x;
subplot(1,2,1), plot(x,y1,x,y2,x,y3)
title('one')
subplot(1,2,2), plot(x,y4,x,y5,x,y6)
title('two')