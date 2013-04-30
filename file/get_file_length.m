% Get the size of a file in lines
%
% The input: 
% file_name - name of input file
% 
% The output: 
% num_lines - number of lines in input file
% 
function num_lines = get_file_length(file_name)

fid = fopen(file_name, 'r'); % open file
num_lines = 0;
while 1
    line=fgetl(fid);
    if (~ischar(line)) %EOF
        break;
    end
    num_lines=num_lines+1;
    if mod(num_lines,1000)==0
        sprintf('scan line %d', num_lines)
    end
end
fclose(fid);
