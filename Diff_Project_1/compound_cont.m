function out = compound_cont(r,t,y_0)
    y = [];
    for i = 1:(length(t))
       y = [y y_0.*exp(r.*t(i))];
    end
    out = y;
end