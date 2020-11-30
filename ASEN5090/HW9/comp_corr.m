function out = comp_corr(data,delay,theta,n,sampprn)
    out = sum(data(delay+1:delay+n).*sampprn.*exp(-1i.*theta));
end