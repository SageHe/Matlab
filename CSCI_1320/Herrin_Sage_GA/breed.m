function out = breed(par1, par2)
% % %Method one
% child = [];
% for i = 1:randi(length(par1))
%     child = [child par1(i)];
% end
% for j = (i + 1):length(par2)
%     child = [child par2(j)];
% end
% out = child;

    
% %Method two
vec = [];
child = [];
for i = 1:length(par1)
    a = randi([1,2]);
    vec = [vec a];
end
for j = 1:length(vec)
    if vec(j) == 1
        child = [child par1(j)];
    else if vec(j) == 2
        child = [child par2(j)];
        end
    end
end
out = child;
end