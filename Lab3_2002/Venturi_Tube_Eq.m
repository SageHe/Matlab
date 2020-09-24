function [ V ] = Venturi_Tube_Eq( auxdeltaP, Tatm, Patm, Aratio )
%Takes in data from the wind tunnel and willreturn the velocity of the flow
R = 287.1; %J/kg/K
Vtop = 2.*auxdeltaP.*R.*Tatm;
Vbot = Patm.*(1-((1/Aratio).^2));
V = (Vtop./Vbot).^(1/2);

end

