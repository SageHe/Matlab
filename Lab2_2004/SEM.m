read_data_alt
N = [1 2]';
STD = [0 0]';
Run_mean = [0 mean(Isp(1:2))]';
std_err = [0 0]';    
ISP = [Isp(1) Isp(2)]';
mat = [N STD Run_mean std_err ISP];
for i = 3:numel(Isp)
    mat(i,1) = i;
    temp = Isp(1:i);
    mat(i,2) = std(temp);
    mat(i,3) = mean(temp);
%     N = numel(temp);
    mat(i,4) = std(temp)/sqrt(numel(temp));
    mat(i,5) = Isp(i);
end
for j = 3:length(mat)
    if 1.96*mat(j,4) < .05*mat(j,3)
        fprintf('It takes %d tests to achieve 95 percent confidence in Isp\n',j)
        break
    end
end
for k = 3:length(mat)
    if 2.23*mat(k,4) < .025*mat(k,3)
        fprintf('It takes %d tests to achieve 97.5 percent confidence in Isp\n',k)
        break
    end
end
for l = 3:length(mat)
    if 2.58*mat(l,4) < .01*mat(l,3)
        fprintf('It takes %d tests to achieve 99 percent confidence in Isp\n',l)
        break
    end
end
acc_99_num = (((2.58*mat(33,2))/.01)^2);
acc_99_num = round(acc_99_num) + 33;
fprintf('It takes %d tests to achieve 99 percent confidence in Isp\n',acc_99_num)
plot(mat(:,1),mat(:,4))
title('SEM vs Trial Number')
xlabel('Number of Trials N')
ylabel('Standard Error of the Mean')