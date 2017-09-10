%function table=loadCellFile_str(infile)
function table=loadCellFile_str(infile)
fin=fopen(infile, 'r');

nLines=0;

while 1
    line=fgetl(fin);
    if (~ischar(line)) %EOF
        break;
    end

    nLines=nLines+1;
    if mod(nLines,200)==0
        disp(nLines);
    end

    tabs=find(line==9); %9 is TAB
    starts=[1 tabs+1];
    ends=[tabs-1 length(line)];

    for i=1:length(starts)
        table{nLines, i}=line(starts(i):ends(i));
    end
end

fclose(fin);
