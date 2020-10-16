function [E] = solvekep(M,e)
E = M;
tol = 0.001;
while (M - (E - e*sin(E))) > tol
    E = E - (M - (E - e*sin(E)))/(-(1 - e*cos(E)));
end
end