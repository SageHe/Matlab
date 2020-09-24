function [ V ] = Pitot_Static_EQ( Air_speed_deltaP, Tatm, Patm )
%Calculated the flow velocity from a pitot static tube
R = 286.9; %J/kg/K
Q = 2.*Air_speed_deltaP.*((R.*Tatm)./Patm);
V = Q.^(1/2);

end

