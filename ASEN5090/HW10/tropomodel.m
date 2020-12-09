function tropo = tropomodel(el,zd)
% tropo = (zd)./sind(el);
m = 1/(sqrt(1 - (cosd(el)/1.001)^2));
tropo = m*zd;
end