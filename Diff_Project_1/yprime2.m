function YP = yprime2(y,t,p)
YP = y*(.03 + .015*sqrt(t - 5)) - 12*p;