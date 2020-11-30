function out = comp_corr(data,delay,tvec,sampprn,IF,fd)
    n = numel(tvec);
    theta = 2*pi*(IF+fd)*tvec;
    out = sum(data(delay+1:delay+n).*sampprn.*exp(-1i.*theta));
end