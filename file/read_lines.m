% Read a set of consecutive lines from a .txt file 
%
% The input: 
% file_name - name of file 
% start_line - first line to read
% end_line - last line to read (default: like first line)
% 
% The ouptut: 
% s - struct with all lines (if one line s is a string)
% 
function s = read_lines(file_name, start_line, end_line, output_file, varargin)


if(~exist('end_line', 'var') || isempty(end_line))
	end_line = start_line;
end
fid = fopen(file_name, 'r'); % open file 
num_lines = end_line - start_line + 1;
for i=1:start_line-1
	dummy = fgetl(fid); 
end
s = cell(num_lines,1); 
for i=1:num_lines
	s{i} = fgetl(fid);
end
fclose(fid); 
if(num_lines == 1)
	s = s{1};
end
if(exist('output_file', 'var'))
	savecellfile(s, output_file); 	
end



