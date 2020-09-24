function out = FOS
failure_prob = input('Input percentage failure probability \n');
F = icdf('normal',failure_prob,4.8,.4);
N = (4.8 - F)/.4;
FOS = 1/(1 - (N*(.4/4.8)));
out = FOS;