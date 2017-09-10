%function words_cell = parse_words_from_line(line_str)
% returns the words separated in a cell
% line_str is a string
function words_cell = parse_words_from_line(line_str)

[temp, rem] = strtok(line_str);
%first cound num words
num_words = 0;
while (min(size(temp))~=0)
    num_words = num_words+1;
    [temp, rem] = strtok(rem);
end

words_cell = cell(num_words, 1);
[word, rem] = strtok(line_str);
for i = 1:num_words
    words_cell{i,1} = word;
    [word, rem] = strtok(rem);
end