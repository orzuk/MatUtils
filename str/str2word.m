% Return one word out of a string, seperated by a certain delimiter
%
% Input:
% delimiter - how to seperate the string into words
% str - the input string
% ind - any positive integer not more than the number of words (can be 'end')
%
% Output:
% w - the word
%
function w = str2word(delimiter, str, ind)
if(iscell(str))
    n = length(str);
    w = cell(n,1);
    for i=1:n
        w{i} = str2word(delimiter, str{i}, ind);
    end
else
    a = strsplit(str, delimiter);
    if(isnumeric(ind))
        if(length(a) >= min(ind))
            w = a{ind};
        else
            w = ''; % return an empty string
        end
    else
        eval(['w = a{' ind '};']);
    end
end
