%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Theodore Trozinski
% Maya West
% ASEN 5044 Estimation; Final Project 
% Purpose: Check measurements, selects common ones if meas length does not
% equal
% Created: April 17, 2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y_,ytru,H] = datcheck(meas,y,Hin)
if isempty(meas) == 1 || isempty(y) == 1
    fprintf('Gap In Data \n')
    y_ = [];
    ytru = [];
    H = [];
elseif numel(meas) == 4 && numel(y) == 4
    y_ = meas(1:3); ytru = y(1:3); H = Hin;
elseif numel(meas) == 8 && numel(y) == 8
    y_ = meas(1:3,:); ytru = y(1:3,:); H = Hin;
elseif numel(meas) == 4 && numel(y) == 8
    vis = meas(end);
    index = find(y==vis);
    if index == 4
        y_ = meas(1:3,1); ytru = y(1:3,1); H = Hin;
    elseif index == 8
        y_ = meas(1:3,1); ytru = y(1:3,2); H = Hin;
    end
elseif numel(meas) == 8 && numel(y) == 4
    vis = y(end);
    index = find(meas==vis);
    if index == 4
        y_ = meas(1:3,1); ytru = y(1:3,1);
        H = Hin(1:3,:);
    elseif index == 8
        y_ = meas(1:3,2); ytru = y(1:3,1);
        H = Hin(4:6,:);
    end
else 
    fprintf('Ruh Roh Shaggy \n')
end