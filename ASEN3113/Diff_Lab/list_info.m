function[list_length,list_max,list_min] = list_info(list)
list_length = length(list);
list_max = max(list);
list_min = min(list);

fprintf('The length of the list = %d \n',list_length)
fprintf('The max of the list = %d \n',list_max)
fprintf('The min of the list = %d \n',list_min)

fprintf('Go Buffs! \n')

end