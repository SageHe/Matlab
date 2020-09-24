function dydt = odefunp4(t,y,Phi)


y = reshape(y,size(Phi));
dydt = Phi*y;

dydt = reshape(dydt,16,1);
end



