% Remove all digits from a string
function s = str2nonums(s)
    for i=0:9
        s = strdiff(s, num2str(i));
    end
end

