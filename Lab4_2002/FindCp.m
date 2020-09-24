function [CNA, CAA, CP] = FindCp(input,scanivalve_pos)


c=0.0889; %metes

for y = 1:12
 
    p(y,:) = input(y,7:22); %Reads in scanivalve pressures. Should it just be 7:22? Why is angel of attack included?
    
        for i = 1:16
        Cp(y,i) = (p(y,i)) / input(y,5); %Calculates coefficient of pressure for each scanivalve port
        end 
        
        Y = [Cp(y,8) Cp(y,9)];
        X = [scanivalve_pos(8,1) scanivalve_pos(9,1)];

        [P S] = polyfit(X,Y,1);
        val1 = polyval(P,3.5*.0254);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Y = [Cp(y,11) Cp(y,10)];
        Z = [scanivalve_pos(11,1) scanivalve_pos(10,1)];

        [Z a] = polyfit(X,Y,1);
        val2 = polyval(Z,3.5*.0254);

        output(y) = (val1+val2)/2;
        
        output = output';
        
        CP(y,:) = [Cp(y,1:9) output(y) Cp(y,10:end)];
        
        scanivalve_pos = [0 .14665; .175 .33075; .35 .04018; .7 .476; 1.05 .49; 1.4 .4774; 1.75 .4403; 2.1 .38325; 2.8 .21875; 3.5 0; 2.8 0; 2.1 0; 1.4 0; 1.05 0; 0.7 .0014; .35 .0175; .0175 .03885];
        scanivalve_pos = scanivalve_pos*0.0254;
        
        for z = 1:((length(CP)) - 1)  %Calcualtes axial and normal coefficient           
%           Cn(z,1) = .5*(Cp(z) + Cp(z+1))*((scanivalve_pos(z+1,1)+scanivalve_pos(z,1))/(2*c));
%           Ca(z,1) = .5*(Cp(z) + Cp(z+1))*((scanivalve_pos(z+1,2)+scanivalve_pos(z,2))/(2*c));
%           Cn(z,1) = trapz(Cp,(scanivalve_pos(:,1)/c));  
%           Ca(z,1) = trapz(Cp,(scanivalve_pos(:,2)/c));  
            Cn(z,1) = (CP(z) + CP(z + 1))*((scanivalve_pos(z+1,1) - scanivalve_pos(z,1))/c);
            Ca(z,1) = (CP(z) + CP(z + 1))*((scanivalve_pos(z+1,2) - scanivalve_pos(z,2))/c);
        end
                
                 a = -sum(Cn);
                 CNA(1,y) = a;
                 b = sum(Ca); 
                 CAA(1,y) = b;
                 
end
end