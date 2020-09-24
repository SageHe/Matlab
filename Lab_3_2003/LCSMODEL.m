%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ASEN 2003 - Lab 3: Model Case
%Purpose: 
%   - Import the derivations for an ideal velocity (V_mod)
%   - Export V_mod
%Created: 2/14/2018
%Modified: 2/19/2018
%Creators:
%   - Lucas Zardini
%   - Sage Herrin
%   - Yang Lee
%   - Idam Isnaeni
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v_mod] = LCSMODEL(r, d, l, theta_exp, w_exp)
%% Derivations 
%original derivation
%v_mod = -r.*w_exp.*(sin(theta_exp)-((cos(theta_exp).*(d-r*sin(theta_exp)))./(l.*cos(asin((d-r.*sin(theta_exp))./l)))));
%updated derivation
v_mod = r.*w_exp.*(-1*sin(theta_exp) - cos(theta_exp).*tan(asin((d-r.*sin(theta_exp))/l)));
end