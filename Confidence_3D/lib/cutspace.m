function nospace = cutspace(str)
space = isspace(str);
index = space==0;
nospace = str(index);
