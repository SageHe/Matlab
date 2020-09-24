r_pt = 149.6e6;
r_at = 227.897396e6;
a_t = 188.748692e6;
mu = 1.32712e11;
delta_v1 = sqrt(((2*mu)/(r_pt)) - (mu/a_t)) - sqrt((mu/r_pt));
delta_v2 = sqrt((mu/r_at)) - sqrt(((2*mu)/r_at) - (mu/a_t));
TOF = pi*sqrt((a_t^3)/mu);