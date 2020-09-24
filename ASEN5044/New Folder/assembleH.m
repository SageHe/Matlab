%% Purpose: determine which Hi matrices to use and assemble them

function [H,included,yNL,observed] = assembleH(t,state)
x = state(1);
xdot = state(2);
y = state(3);
ydot = state(4);

RE = 6378;
wE = 2*pi/86400;


included = [];
for i = 1:12
    theta0 = (i-1)*pi/6;
    
    xs = RE*cos(wE*t + theta0);
    ys = RE*sin(wE*t + theta0);
    xdots = -wE*RE*sin(wE*t + theta0);
    ydots = wE*RE*cos(wE*t + theta0);
    
    theta = atan2(ys,xs);
    phi = atan2((y-ys),(x-xs));
    
    if theta >= pi/2
        if phi <= theta-3*pi/2
            included = [included i];
        end
    end
    if theta <= -pi/2
        if phi >= theta+3*pi/2 
            included = [included i];
        end
    end
    
    if phi >= -pi/2+theta && phi <= pi/2 + theta
        included = [included i];
        phi;
    end

    h1x1 = (x-xs)/sqrt((x-xs)^2 + (y-ys)^2);
    h2x1 = (xdot-xdots)/sqrt((x-xs)^2 + (y-ys)^2) - ...
        ((x-xs)*(xdot-xdots) + (y-ys)*(ydot-ydots))*(x-xs)/(((x-xs)^2 + (y-ys)^2)^(3/2));
    h3x1 = (y-ys)/((y-ys)^2 + (x-xs)^2);

    h1x2 = 0;
    h2x2 = (x-xs)/sqrt((x-xs)^2 + (y-ys)^2);
    h3x2 = 0;
    
    h1x3 = (y-ys)/sqrt((x-xs)^2 + (y-ys)^2);
    h2x3 = (ydot-ydots)/sqrt((x-xs)^2 + (y-ys)^2) - ...
        ((x-xs)*(xdot-xdots) + (y-ys)*(ydot-ydots))*(y-ys)/(((x-xs)^2 + (y-ys)^2)^(3/2));
    h3x3 = (x-xs)/((y-ys)^2 + (x-xs)^2);
    
    h1x4 = 0;
    h2x4 = (y-ys)/sqrt((x-xs)^2 + (y-ys)^2);
    h3x4 = 0;
    
    Hi{i} = [[h1x1 h2x1 h3x1]',[h1x2 h2x2 h3x2]',[h1x3 h2x3 h3x3]',[h1x4 h2x4 h3x4]'];
    
    % Nonlinear measurements
    rho = sqrt((x-xs)^2 + (y-ys)^2);
    rhodot = ((x-xs)*(xdot-xdots) + (y-ys)*(ydot-ydots)) / rho;
    yNL{i} = [rho; rhodot; phi];
    
end
yNL = cat(2,yNL{included});
% Assemble H matrices
H = cat(1,Hi{included});

observed = cat(1,yNL,cat(2,included));

end