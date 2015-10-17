%function table=loadcellfile(infile)
% Faster version
function table=loadcellfile(infile, to_mat, delimiter, num_lines)

if(~exist('to_mat', 'var') || isempty(to_mat))
    to_mat = 1; 
end
%if(~exist('delimiter', 'var') || isempty(delimiter))
%    delimiter = 9; 
%end


%find number of columns
fin=fopen(infile, 'r');
%read a couple of lines to try and avoid comment lines etc.
num_of_cols=0;
for nLines=1:20
    line=fgetl(fin);
    if (~ischar(line)) %EOF
        break;
    end
    num_of_cols=max(num_of_cols,length(find(line==9))+1);
end
fclose(fin);

%find number of rows
fin=fopen(infile, 'r');
chunksize = 1e6; % read chuncks of 1MB at a time
nLines = 0;
while ~feof(fin)
    ch = fread(fin, chunksize, '*uchar');
    if isempty(ch)
        break
    end
    nLines = nLines + sum(ch == sprintf('\n'));
end
fclose(fin); 

table=cell(nLines, num_of_cols); cur_num_of_cols = num_of_cols;
fin=fopen(infile, 'r');
nLines=0;
while 1
    line=fgetl(fin);
    if (~ischar(line)) %EOF
        break;
    end

    nLines=nLines+1;
    if mod(nLines,1000)==0
        disp(['read ',num2str(nLines/1000),'K lines']);
    end
    
    if(~exist('delimiter', 'var') || isempty(delimiter))
        tabs=find(line==9 | line==' '); %9 is TAB. Look for tabs or spaces 
    else
        tabs = find(line == delimiter);
    end
    starts=[1 tabs+1];
    ends=[tabs-1 length(line)];
    if(length(starts) > cur_num_of_cols)
        sprintf('adding %d columns at %d row', length(starts), nLines)
        cur_num_of_cols = length(starts);
    end

    for i=1:length(starts)
        dum=sscanf(line(starts(i):ends(i)),'%f');
        if isempty(dum)
            table{nLines, i}=line(starts(i):ends(i));
        else
            if(to_mat)
                dum=str2num(line(starts(i):ends(i)));
                if isempty(dum)
                    table{nLines, i}=line(starts(i):ends(i));
                else
                    table{nLines, i}=dum;
                end
            else
                table{nLines, i}=line(starts(i):ends(i));
            end
        
        end
    end
end

fclose(fin);
