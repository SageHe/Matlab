close 

figure
plot(ell_50(1,:)+BdotT150,ell_50(2,:)+BdotR150,'.')
hold on
plot(BdotT150,BdotR150,'x','MarkerSize',30)
plot(ell100(1,:)+BdotT1100,ell100(2,:)+BdotR1100,'.')
plot(BdotT1100,BdotR1100,'x','MarkerSize',30)
plot(ell150(1,:)+BdotT1150,ell150(2,:)+BdotR1150,'.')
plot(BdotT1150,BdotR1150,'x','MarkerSize',30)
plot(ell200(1,:)+BdotT1200,ell200(2,:)+BdotR1200,'.')
plot(BdotT1200,BdotR1200,'x','MarkerSize',30)
grid on
grid minor
xlabel('$\hat{T}$ Component','Interpreter','latex')
ylabel('$\hat{R}$ Component','Interpreter','latex')
title('B-Plane Estimate with $3\sigma$ Covariance Ellipses','Interpreter','latex')
legend('50 day $3\sigma$ covariance ellipse','50 day Estimate','100 day $3\sigma$ covariance ellipse','100 day estimate','150 day $3\sigma$ covariance ellipse','150 day estimate','200 day $3\sigma$ covariance ellipse','200 day estimate','Interpreter','latex')

set(gca,'FontSize',20)
