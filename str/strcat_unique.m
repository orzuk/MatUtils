% Concatenate two strings, but remove duplicate substrings
function s = strcat_unique(s1, s2)

s = strcat(s1, strdiff(s2,s1));
