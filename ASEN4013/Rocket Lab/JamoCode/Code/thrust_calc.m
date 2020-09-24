function [ThSM_en, Cs, m_dot] = thrust_calc(Pa, Pc, Ae, rho_s, burn_rate, A_burn, AR_sup, delta_Vol) % - used eqn 10.57 for m_dot
    Pc_en = Pc * 145.038; % [psi] chamber pressure

    %% CEA RUN 
    ERROR = 0;
    try % tests if there is any CEA output
        % OUTPUT1 GIVES VALUES IN THE CHAMBER
        % OUTPUT3 GIVES VALUES AT NOZZLE EXIT
        [Output1, Output3] = RUN_CEA(Pc_en, AR_sup);
    catch
       ERROR = 1;
    end
    
    if ERROR == 1 % sets Mach and alpha to zero if output DNE
        alpha = 0;
        Mach = 0;
        Pe = Pa;
        Cs = 0;
        rho_g = 0;
    else
        alpha = Output3.a; % [m/s] sonic velocity
        Mach = Output3.Mach; % Mach number
        Pe = Output3.P * 1e5; % [Pa] nozzle exit pressure
        Cs = Output3.Cstar; % [m/s] characteristic velocity
        rho_g = Output1.rho; % [kg/m^3] propellant gas density, chamber
    end

    %% MASS FLOW CALCULATION
    %m_dot =(rho_s - rho_g)*burn_rate*A_burn; % [kg/s] propellant mass flow rate
    m_dot = delta_Vol*rho_s;

    %% THRUST CALCULATION
    Ve = alpha*Mach;
    ThSM = (m_dot*Ve) + (Pe - Pa)*Ae; % [N] thrust SI
    ThSM_en = ThSM * 0.224809; % [lbf] imperial thrust to match curve data
