%% 4, same process done for each part of 4
DCM = [0 -1 0;-1 0 0;0 0 -1]*[0 1/sqrt(2) -1/sqrt(2);1 0 0;0 -1/sqrt(2) -1/sqrt(2)];
[V,D] = eig(DCM);
axis = V(:,3); % different column of V is chosen depending on column that corresponds to same column of D that contains a 1

