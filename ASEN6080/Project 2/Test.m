opts = odeset('RelTol',1e-12,'AbsTol',1e-12);
load('Project2_Prob2_truth_traj_50days.mat')

tspan = Tt_50;
state = Xt_50(1,1:7);

[time,X] = ode45(@(t,state) TBG_SRP_ode(t,state),tspan,state,opts);
