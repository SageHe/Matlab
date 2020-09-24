function [sigma_BR,Bw_BR] = att_err(sigma_BN,Bw_BN,RN,Nw_RN)
BN = MRP2C(sigma_BN);
Bw_RN = BN*Nw_RN;
Bw_BR = Bw_BN - Bw_RN;
sigma_RN = C2MRP(RN);
sigma_BR = addMRP(-sigma_RN,sigma_BN);
if norm(sigma_BR) > 1
    sigma_BR = -sigma_BR/(sigma_BR'*sigma_BR);
end
end
 