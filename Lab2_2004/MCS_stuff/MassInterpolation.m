% Split the mass array up based on time separation between datapoints. 
% Use interpolation between consecutive mass points to account for spots
% in between the time samples.
% Having thrust dependent on time, or thrust as a function of time, will
% allow ODE45 to use it. 

function m = MassInterpolation(t,m,mtimes,mfinal)

% check to see if t > the largest value of time the mass data spans. If so,
% then the mass = the final mass of the rocket.
if mtimes(end-1) < t
    m = mfinal;

else
    

% Loop through to find the index of mtimes that corresponds to t.
j = 1;
while (mtimes(j) < t)
    j = j+1;
end

% interpolate to find mass
%m = abs(((m(j+1) - m(j))/((mtimes(j+1) - mtimes(j)))) * (t - mtimes(j)) + m(j));
minterp = (m(j+1)+m(j))/2;

if minterp > mfinal
    m = minterp;
else
    m = mfinal;
end

end

end