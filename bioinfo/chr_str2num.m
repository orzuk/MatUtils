% Convert string to number with special attention to x,y and m chroms
% Note that this depend on the organism. Currently this is good for mammals
% Input: 
% chr_str - string with the chromosome name 
% organism_str - organism (optional. Default is human)
% Output: 
% chr_num - number of chromosome 
%
function chr_num = chr_str2num(chr_str, organism_str, varargin)

if(iscell(chr_str))
    chr_num = zeros(length(chr_str),1);
    for i=1:length(chr_str)
        if(nargin == 1) % no organism means human
            chr_num(i) = chr_str2num(chr_str{i});
        else
            chr_num(i) = chr_str2num(chr_str{i}, organism_str);
        end
    end
    return;
end

if(~isempty(strfind(chr_str, 'random'))) % New! random is considered as an error
    chr_num = -1; return;
end

if(my_width(chr_str) == 1) % only one to convert
    chr_str = strdiff(chr_str, 'chr'); % remove possible prefix
end
if(nargin == 1) % no organism means human
    if(strcmp(chr_str, 'x') || strcmp(chr_str, 'X'))
        chr_num = 23; return;
    end
    if(strcmp(chr_str, 'y') || strcmp(chr_str, 'Y'))
        chr_num = 24; return;
    end
    if(strcmp(chr_str, 'm') || strcmp(chr_str, 'M'))
        chr_num = 25; return;
    end
    if(strcmp(chr_str, '-')) % New! used to handle errors
        chr_num = -1; return;
    end
    chr_num = str2num(chr_str);
else
    Assign24MammalsBasicGlobalConstants;
    organism_str = genome_ver_to_organism(organism_str); % enable also input as genome version
    organism_ind = -1;
    for i_org_str=1:length(genome_versions_vec)
        if(strcmp(organism_str, genome_versions_vec{i_org_str}{1}))
            organism_ind = i_org_str;
            break;
        end
    end
    organism_chr = genome_versions_vec{organism_ind}{3};
    %    ind = find(strcmp(upper(organism_str), MANUEL_SPECIES_ORDER));

    switch chr_str
        case {'x', 'X'}
            chr_num = organism_chr; return;
        case {'y', 'Y'}
            chr_num = organism_chr+1; return;
        case {'m', 'M'}
            chr_num = organism_chr+2; return;
    end
    chr_num = str2num(chr_str);
end

if(isempty(chr_num))
	chr_num = -1; 
end


