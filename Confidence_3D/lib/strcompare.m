function out = strcompare(str1,str2)

% str1: string matrix to be compared 
% str2 is a string vector. 
% scan through the rows of str1 to find if any rows equal to str2 (space is...
% not considered to compare)
% if there is some rows equal, return out=1


out = 0;%defaut
length1 = size(str1,1);

for it = 1:length1
    temp = cutspace(str1(it,:));
    if length(temp)==length(str2)
        out = strcmp(temp,str2);
        if out == 1
            break;
        end;
    end
end




