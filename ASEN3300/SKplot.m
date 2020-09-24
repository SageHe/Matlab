function[w0,Q] = SKplot(R,C1,C2)
w0 = 1/(R*sqrt(C1*C2));
Q = .5*sqrt(C1/C2);
H = tf([w0^2],[1 (w0/Q) w0^2]);
opts = bodeoptions;
opts.FreqScale = 'linear';
opts.FreqUnits = 'Hz';
w = [0:60000];
bode(H,w,opts);
title('Phase and Magnitude vs Frequency')
end