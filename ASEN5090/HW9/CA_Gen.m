function [out] = CA_Gen(in,length)
G1 = ones(1,10);
G2 = ones(1,10);
% if in == 19
%     PS = xor(G2(3), G2(6));
% elseif in == 25
%     PS = xor(G2(5),G2(7));
% elseif in == 5
%     PS = xor(G2(5),G2(7));
% end
for i = 1:length
    temp = [G2(2) G2(3) G2(6) G2(8) G2(9) G2(10)];
    G2out(i) = G2(end);
if in == 19
    PS = xor(G2(3), G2(6));
elseif in == 25
    PS = xor(G2(5),G2(7));
elseif in == 5
    PS = xor(G2(1),G2(9));
elseif in == 31
    PS = xor(G2(3),G2(8));
elseif in == 1
    PS = xor(G2(2),G2(6));
elseif in == 22
    PS = xor(G2(6),G2(9));
elseif in == 32
    PS = xor(G2(4),G2(9));
elseif in == 10
    PS = xor(G2(2),G2(3));
end
%     PS = xor(G2(3),G2(6));
    PRN(i) = xor(PS,G1(end));
    G1out(i) = G1(end);
    G1 = [xor(G1(3),G1(10)) G1(1:9)];
    G2 = [rem(sum(temp),2) G2(1:9)];
end
out = PRN;