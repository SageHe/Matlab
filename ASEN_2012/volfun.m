function yprime = odefun(t,y_0)
global g C_d rho_air_amb Vol_bottle p_amb gamma rho_water D_throat D_bottle R MBottle C_D p_gage Vol_water_i T_air_i v0 theta x0 y0 l_s tspan

yprime = 0;
yprime = C_d*(pi*(D_throat./2)^2)*sqrt((2./rho_water)*((p_gage + p_amb)*(.001./y_0)^gamma - p_amb));

end

