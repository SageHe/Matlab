function[dydt] = HW5ode(t,y)
global g m rad I_x I_y I_z k alpha eta f f1 f2 f3 f4 X Y Z K1_roll K2_roll K1_pitch K2_pitch K_y feedback
% Monday 11:25, 5 deg. at 1.313 or wed 12:53, 5 degrees at .5 seconds
    u_E = y(1);
    v_E = y(2);
    w_E = y(3);
    aero_forces = -eta*norm([u_E v_E w_E])*[u_E v_E w_E];
    X = aero_forces(1);
    Y = aero_forces(2);
    Z = f + aero_forces(3);
    p = y(4);
    q = y(5);
    r = y(6);
    psi = y(7);
    theta = y(8);
    phi = y(9);
    N = y(10);
    E = y(11);
    D = y(12);
    phi_d = y(13);
    theta_d = y(14);
    if t(end) < 1.313
%         L_c = (rad/sqrt(2)*(f1 + f2 - f3 - f4));
%         M_c = (rad/sqrt(2))*(f1 + f4 - f2 - f3);
%         N_c = (k*(f2 + f4 - f1 - f3));
        L_c = -K1_roll*p - K2_roll*(phi - phi_d);
        M_c = -K1_pitch*q - K2_pitch*(theta - theta_d);
        N_c = -K_y*r;

    else
        L_c = -K1_roll*p - K2_roll*phi;
        M_c = -K1_pitch*q - K2_pitch*theta;
        N_c = -K_y*r;
    end
    G = [[L_c M_c N_c] + -alpha*norm([p q r])*[p q r]];
    L = G(1);
    M = G(2);
    N = G(3);
        

    dudt = X/m - g*sin(theta) - q*w_E + r*v_E;
    dvdt = Y/m + g*cos(theta)*sin(phi) - r*u_E + p*w_E;
    dwdt = Z/m + g*cos(phi)*cos(theta) - p*v_E + q*u_E;
    dpdt = (L - q*r*(I_z - I_y))*(1/I_x);
    dqdt = (M - r*p*(I_x - I_z))*(1/I_y);
    drdt = (N - p*q*(I_y - I_x))*(1/I_z);
    dphidt = p + (q*sin(phi) + r*cos(phi))*tan(theta);
    dthetadt = q*cos(phi) - r*sin(phi);
    dpsidt = (q*sin(phi) + r*cos(phi))*sec(theta);

    dydt(1) = dudt;
    dydt(2) = dvdt;
    dydt(3) = dwdt;
    dydt(4) = dpdt;
    dydt(5) = dqdt;
    dydt(6) = drdt;
    dydt(7) = dpsidt;
    dydt(8) = dthetadt;
    dydt(9) = dphidt;
    dydt(10) = u_E*cos(theta)*cos(psi) + v_E*(sin(phi)*sin(theta)*cos(psi) - cos(phi)*sin(psi)) + w_E*(cos(phi)*sin(theta)*cos(psi) + sin(phi)*sin(psi));
    dydt(11) = u_E*cos(theta)*sin(psi) + v_E*(sin(phi)*sin(theta)*sin(psi) + cos(phi)*cos(psi)) + w_E*(cos(phi)*sin(theta)*sin(psi) - sin(phi)*cos(psi));
    dydt(12) = -u_E*sin(theta) + v_E*sin(phi)*cos(theta) + w_E*cos(phi)*cos(theta);
    dydt(13) = 0;
    dydt(14) = 0;

    dydt = dydt';

% end

