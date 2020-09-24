function out = isPalindrome(string,n)
string = lower(string);
new_string = [];
for i = 1:length(string)
    new_string(i) = string(end - (i - 1));
end
new_string = char(new_string);
if isequal(new_string, string) == 1
    new_string(1:n) = [];
    for j = 1:n
        new_string(end) = [];
    end
else
    out = isequal(new_string, string);
    return;
end
out = new_string;
%out = isequal(new_string, string);
end
 
