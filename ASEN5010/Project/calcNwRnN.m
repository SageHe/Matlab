function [NwRnN] = calcNwRnN(t)
global thetadotLMO
HN = calcHN(t);
w_hn = [0 0 thetadotLMO]';
NwRnN = HN'*w_hn;
end


