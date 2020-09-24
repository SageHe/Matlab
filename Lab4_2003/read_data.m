%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ASEN 2003 Lab 4: BALANCED AND UNBALANCED WHEEL
%
% Created:  03/05/2018 - Charles Puskar
% Modified: 03/05/2018 - Charles Puskar
%
% read_data.m
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [time, theta_out, w_out] = read_data(filename)
    % opens document, reads first two commented lines
    fid = fopen(filename);
    fgetl(fid);
    fgetl(fid);
    
    % reads data into matrix
    A = fscanf(fid,'%f',[3, Inf]);

    % closes file
    fclose(fid);
    
    % separates data into arrays of components
    time = A(1,:)';
    theta = A(2,:)';
    w = A(3,:)';
    
    remove = [];
    for i = 1:length(theta)
        if theta(i)<-1 || theta(i)>15
            remove = [remove,i];
        end
    end
    
    remove
    
    time(remove) = [];
    theta(remove) = [];
    w(remove) = [];
    
    theta_out = theta;
    w_out = w;
end