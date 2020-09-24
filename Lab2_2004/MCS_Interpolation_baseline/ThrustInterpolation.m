% Split the Thrust array up based on time separation between datapoints so
% that it's dependent on time. 
% Use interpolation between consecutive Thrust points to account for spots
% in between the time samples.
% Having thrust dependent on time, or thrust as a function of time, will
% allow ODE45 to use it. 

function T = ThrustInterpolation(t,F)

% set the sampling rate
% Multiply the sampling rate by an experimentally determined constant so that
% the thrust runs out when the mass runs out.

samplingrate = 0.001652*0.55;

% make an array of time values that correspond to each point in the experimental thrust
% array
Ftimes = zeros(1,length(F));
for i = 1:length(F)
    Ftimes(i) = i*samplingrate;
end


% now determine the Thrust, T. If t > the time that the experimental data
% spans, then T = 0. If not, the T = some interpolated value from the
% thrust data corresponding to t.

% use end-1 because if only use end, then the last value will try to index
% end+1 in the force array, which doesn't exist.
if Ftimes(end-1) < t
    T = 0;
    
else
    % Find the index in Ftimes that corresponds to t
    j = 1;
    % check if the given time t is within the range of the time that the thrust
    % experiment took place
    while (Ftimes(j) < t)
        j = j+1;
    end
    
    % now eliminate bad experimental data
    if F(j)<0 %|| F(j)>200
        T = 0;

    % interpolate to find Thrust
    else        
        %T = abs(((F(j+1) - F(j))/((Ftimes(j+1) - Ftimes(j)))) * (t - Ftimes(j)) + F(j));
        T = (F(j+1) + F(j))/2;
    end
end

end