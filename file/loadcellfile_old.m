% Load a tab-delimited .txt file into a cel array (Written by Assif Yitzhaky)
%
% Input:
% infile - input file in .txt format (tab-delimited)
% to_mat - flag saying if to convert strings to mat (default or not)
%
% Output:
% table - a cell array containing the .txt file
%
function table=loadcellfile(infile, to_mat, varargin)

fin=fopen(infile, 'r');
table = {}; % start with an empty cell
nLines=0;
if(~exist('to_mat', 'var') || isempty(to_mat))
    to_mat = 1; 
end

while 1
    line=fgetl(fin);
    if (~ischar(line)) %EOF
        break;
    end    
    nLines=nLines+1;
    if(nLines == 309) 
        ccc = 2312
    end
    if mod(nLines,500)==0
        sprintf('load line %d', nLines)
    end
    
    tabs=find(line==9); %9 is TAB
    starts=[1 tabs+1];
    ends=[tabs-1 length(line)];
    
    for i=1:length(starts)
        dum=sscanf(line(starts(i):ends(i)),'%f');
        table{nLines, i}=line(starts(i):ends(i));
        if(~isempty(dum) && to_mat)
            dum=str2num(line(starts(i):ends(i)));
            if ~isempty(dum)
                table{nLines, i}=dum;
            end
        end
    end
end

fclose(fin);
