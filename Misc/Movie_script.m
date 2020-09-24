clear all;close all;clc
fid = fopen('Movie_list.txt','r');
text = textscan(fid,'%s','Delimiter','\n');
text = string(text{1});
fclose(fid);

ans = 'no';
while ans == 'no'
num = randi(length(text));
movie = text(num);
while movie == 'blank'
  num = randi(length(text));
  movie = text(num);
end
if num < 23
    fprintf('From Kaitlyns Dope Ass Picks, Watch %s\n',movie)
else
    fprintf('From Sages Fucking Choice Fucking Cinema Fucking Selections, Watch %s\n',movie)
end

ans = input('Do you want to watch this?\n');
ans = lower(ans);
if strcmp(ans,'yes')
    break
end
end

text(num) = 'blank';

newfid = fopen('Movie_list.txt','w');
fprintf(newfid,'%s\n',text);
fclose(newfid);

watched_fid = fopen('Watched_Movies.txt','a');
fprintf(watched_fid,'%s was watched on %s\n',movie,datestr(now));
fclose(watched_fid);

fprintf('Nice Choice\n')


    
