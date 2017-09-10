% Convert region type numbering coding to string
function region = region_type_str2num(region_str)

Assign24MammalsGlobalConstants;

if(~iscell(region_str)) % one string
    region = strmatch(region_str, genome_types);
else  % many strings in a cell array 
    num_regions = length(region_str); % get number of regions 
    region = zeros(num_regions, 1);
    for i=1:length(genome_types)
        region(strmatch(genome_types{i}, region_str)) = i;
    end
end
% genome_types = {'intergenic', 'exon', 'utr3', 'intron', 'utr5', 'promoter'}; % matching order to python file ReadHg17File.py
% INTERGENIC=1; EXON=2; UTR3=3; INTRON=4; UTR5=5; PROMOTER=6;  % numeric values to different genome regions

