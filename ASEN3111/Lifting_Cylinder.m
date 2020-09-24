%
% This function plots the streamlines for the flow around a cylinder with
% circulation.  The flow is constructed using a uniform flow, a dipole
% (source and sink at the same location), and a vortex.  Output is a plot
% showing the streamlines around the cylinder.

clear all;clc;clf;close all

%% Define Domain
xmin=-3.5;
xmax=1.5;
ymin=-3.0;
ymax=1.5;

%% Define Number of Grid Points
nx=100; % steps in the x direction
ny=100; % steps in the y direction

%% Create mesh over domain using number of grid points specified
[x,y]=meshgrid(linspace(xmin,xmax,nx),linspace(ymin,ymax,ny));

%% Define Flow Parameters

U0 = 5.0;       % Uniform flow velocity
Gamma = 8*pi;   % 8*pi vortex strength
xGamma = -1.0;  % location of vortex
yGamma = -1.0;
K = 5.0;        % dipole strength
xK = -1.0;      % location of dipole
yK = -1.0;

%% Define a function which calculates the radius.
% Center of circle = (x1,y1)
radius=inline('sqrt((x-x1).^2+(y-y1).^2)','x','y','x1','y1');

%% Calculate psi for uniform stream (Eq. 3.55; pg. 310)
Psi_U0 = U0*y;

%% Calculate psi for doublet or dipole (Eq. 3.87; pg. 310)
Psi_K = K*sin(atan2(y-yK,x-xK))./(radius(x,y,xK,yK));

%% Calculate psi for vortex (Eq. 3.114; pg. 310)
Psi_Gamma = Gamma/(2*pi)*log(radius(x,y,xGamma,yGamma));

%% Add all streamfunctions together
StreamFunction = Psi_U0 - Psi_Gamma - Psi_K;

%% Determine color levels for contours
levmin = StreamFunction(1,nx); % defines the color levels -> trial and error to find a good representation
levmax = StreamFunction(ny,nx/2);
levels = linspace(levmin,levmax,50)';

%% Plot streamfunction at levels
contour(x,y,StreamFunction,levels)

%% Plot cylinder on same graph
theta=linspace(0,2*pi);

figure(1)
hold on;
plot(xGamma+cos(theta),yGamma+sin(theta),'k');
hold off;

%% Adjust axis and label figure
axis equal
axis([xmin xmax ymin ymax])
ylabel('y')
xlabel('x')
