list = [];
num = [];
while length(num) ~= 4
    num = input('Input the desired number of the list or type done in single quotes \n');
    if length(num) ~= 4
        list = [list num];
    end
end
clear num
[list_length,list_max,list_min] = list_info(list);