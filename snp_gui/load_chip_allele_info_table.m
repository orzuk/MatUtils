function load_chip_allele_info_table(chip_name, genome_build, annot_file)

raw_DB=fullfile('..','raw_database');

fin=fopen(fullfile(raw_DB, annot_file), 'r');

titles={'Probe Set ID' 'dbSNP RS ID' 'Chromosome' 'Physical Position' 'Allele A' 'Allele B' 'Strand'};

line1=fgetl(fin);
tabs=strfind(line1,'","');
starts=[2 tabs+3];
ends=[tabs-1 length(line1)-1];

for i=1:length(starts)
    dum=sscanf(line1(starts(i):ends(i)),'%f');
    if isempty(dum)
        line_cell{i}=line1(starts(i):ends(i));
    else
        dum=str2num(line1(starts(i):ends(i)));
        if isempty(dum)
            line_cell{i}=line1(starts(i):ends(i));
        else
            line_cell{i}=dum;
        end
    end
end

snp_table=cell(1,length(titles));
snp_table(1,:)=titles;

cols=zeros(1,length(titles));
for i=1:length(cols)
    cols(i)=find(strcmpi(titles{i},line_cell));
end

nLines=1;

while 1
    line1=fgetl(fin);
    if (~ischar(line1)) %EOF
        break;
    end

    nLines=nLines+1;
    if mod(nLines,200)==0
        disp(nLines);
    end
    
    tabs=strfind(line1,'","');
    starts=[2 tabs+3];
    ends=[tabs-1 length(line1)-1];

    for i=1:length(starts)
        dum=sscanf(line1(starts(i):ends(i)),'%f');
        if isempty(dum)
            line_cell{i}=line1(starts(i):ends(i));
        else
            dum=str2num(line1(starts(i):ends(i)));
            if isempty(dum)
                line_cell{i}=line1(starts(i):ends(i));
            else
                line_cell{i}=dum;
            end
        end
    end
    snp_table(nLines, :)=line_cell(cols);
end

fclose(fin);

save(fullfile(raw_DB,[chip_name '_allele_info_table_' genome_build '.mat']),'snp_table');

