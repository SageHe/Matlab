%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ASEN 2003 - Lab 3: Data Import
%Purpose: 
%   - Import data
%   - Export three variables (theta, w, V)
%Created: 2/14/2018
%Modified: 2/26/2018
%Creators:
%   - Lucas Zardini
%   - Sage Herrin
%   - Yang Lee
%   - Idam Isnaeni
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [theta_exp,w_exp,v_exp] = LCSDATA(fname)
%% Import data
%import data into variable 'Data'
Data = importdata(fname);
%covert 'Data' from struct to a cell
Data = struct2cell(Data);

%exploit data values
%nums - a cell filled with numbers
nums = Data{1};
%time - the time column
time = nums(:,1);
%theta_Exp - the theta (in degrees) column
theta_exp = nums(:,2); %degrees

%%%
%zeros wheel position 
%find a wheel position between 0 and 360 for a starting point
i = theta_exp(1);
j=0;
while i < 0 || i > 360
    i = i-360;
    j=j+1;
end
theta_exp = theta_exp - j*360;
%%%

%conversion from degree to rad
theta_exp = ((theta_exp*pi)/180); % convert to rad
%slide position
slide_position = nums(:,3); %mm
%w_exp - angular velocity
w_exp = nums(:,4); %deg/s
%converts to rad/s
w_exp = ((w_exp*pi)/180); %rad/s
%v_exp - experimental velocity
v_exp = nums(:,5); %mm/s


end