function [RnN] = calcRnN(t)
RH = [-1 0 0;0 1 0;0 0 -1];
HN = calcHN(t);
RnN = RH*HN;
end