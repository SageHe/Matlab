function out = calcECEF2ENU(phi,lambda)
out = [-sind(lambda) cosd(lambda) 0;...
            -sind(phi)*cosd(lambda) -sind(phi)*sind(lambda) cosd(phi);...
            cosd(phi)*cosd(lambda) cosd(phi)*sind(lambda) sind(phi)];
end