% Convert region type numbering coding to string
function region_str = region_type_num2str(region)

Assign24MammalsGlobalConstants;

num_regions = length(region); % get number of regions 
if(num_regions == 1)
    region_str = genome_types{region};
else
    region_str = cell(num_regions, 1);
    for i=1:num_regions
        region_str{i} = genome_types{region(i)};
    end
end

% genome_types = {'intergenic', 'exon', 'utr3', 'intron', 'utr5', 'promoter'}; % matching order to python file ReadHg17File.py
% INTERGENIC=1; EXON=2; UTR3=3; INTRON=4; UTR5=5; PROMOTER=6;  % numeric values to different genome regions
