function [BN,BdotR1,BdotR2,BdotT1,BdotT2,LTOF] = Bplane_calcs(rinf,vinf)
    mu_E = 398600.432896939;
    ecc = (1/mu_E)*((norm(vinf)^2 - (mu_E/norm(rinf)))*rinf - dot(rinf,vinf)*vinf);
    Phat = ecc/norm(ecc);
    h = cross(rinf,vinf);
    What = h/norm(h);
    Qhat = cross(What,Phat);
    
    p = (norm(h)^2)/mu_E;
    
%     a = p/(1 - norm(ecc));
    a = -mu_E/(2*((norm(vinf)^2/2) - (mu_E/norm(rinf))));
    
    b = abs(a)*sqrt(norm(ecc)^2 - 1);
    
    Shat = vinf/norm(vinf);
    Nhat = [0 0 1]';
    That = (cross(Shat,Nhat))/norm(cross(Shat,Nhat));
    
    Rhat = cross(Shat,That);
    
    BN = [Shat;That;Rhat];
    
    BdotR1 = dot(rinf,Rhat);
    BdotT1 = dot(rinf,That);
    
    B = b*(cross(Shat,What));
    
    BdotR2 = dot(B,Rhat);
    BdotT2 = dot(B,That);
    
    cosnu = (rinf/norm(rinf))*Phat';
    
    num = a*(1 - norm(ecc)^2);
    denom = (1 + norm(ecc)*cosnu);
    
    f = acosh(1 + (norm(vinf)^2/mu_E)*(num/denom));
    
    LTOF = (mu_E/norm(vinf)^3)*(sinh(f) - f);
end