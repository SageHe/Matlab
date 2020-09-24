%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2003 Lab 4: BALANCED AND UNBALANCED WHEEL
%
% Created:  02/28/2018 - Charles Puskar
% Modified: 02/28/2018 - Charles Puskar
%
% model_2.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = model_2(theta)

    % defines constants of apparatus
    m_cyl = 11.7; % kg
    m_trail = 0.7; %kg
    r_cyl = 0.235; % m
    k_cyl = 0.203; % m
    I_cyl = m_cyl*k_cyl^2; % kg*m^2
    beta = 5.5*pi/180; % rad
    g = 9.81; % m/s^2
    M = -0.7; % N*m
    
    w = sqrt(((m_cyl+m_trail).*g.*r_cyl.*theta.*sin(beta)+M.*theta)./(0.5.*(m_cyl+m_trail).*r_cyl.^2 + 0.5.*I_cyl));
    
end