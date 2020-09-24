function [w_RN] = calc_wRN(sigma_RN,sigmadot_RN)
w_RN = (4/(1+sigma_RN'*sigma_RN)^2)*((1-sigma_RN'*sigma_RN)*eye(3) - 2*tilde(sigma_RN) + 2*sigma_RN*sigma_RN')*sigmadot_RN;
end