%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2003 Lab 4: BALANCED AND UNBALANCED WHEEL
%
% Created:  02/28/2018 - Charles Puskar
% Modified: 02/28/2018 - Charles Puskar
%
% model_4.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function w = model_4(theta)

    % defines constants of apparatus
    m_cyl = 11.7; % kg
    m_trail = 0.7; %kg
    m_extra = 3.4; % kg
    r_cyl = 0.235; % m
    k_cyl = 0.203; % m
    I_cyl = m_cyl*k_cyl^2; % kg*m^2
    beta = 5.5*pi/180; % rad
    r_to_extra = 0.178; % m
    r_extra = 0.019; % m
    g = 9.81; % m/s^2
    I_extra = 0.5*m_extra*r_extra^2; % kg*m^2
    M = -0.85; % N*m
    
    w = sqrt(((m_cyl+m_trail).*g.*r_cyl.*theta.*sin(beta)+M.*theta+m_extra.*g.*(r_to_extra.*theta.*sin(beta)+r_to_extra-r_to_extra.*cos(theta+beta)))./(0.5.*(m_cyl+m_trail).*r_cyl.^2 + 0.5.*m_extra.*(((r_cyl+r_to_extra.*cos(theta)).^2)+(r_to_extra.*sin(theta)).^2) + 0.5.*I_cyl + 0.5.*I_extra));
    
end