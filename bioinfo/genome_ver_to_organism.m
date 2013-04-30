% A simple function telling us for a given genome version what is the
% organism (currently works only for human and mouse)
%
% Input:
% genome_version - a string containing the (EXACT!) genome version as in ucsc
% Output:
% organism_str - string with the organism's name (in uppercase)
%
function organism_str = genome_ver_to_organism(genome_version)

Assign24MammalsBasicGlobalConstants; % read the vector of genome versions

if(iscell(genome_version))
    num_org = length(genome_version);
    organism_str_vec =cell(num_org,1);
    for i=1:num_org
        organism_str_vec{i} = genome_ver_to_organism(genome_version{i});
    end
    organism_str = organism_str_vec;
else
    for i=1:length(genome_versions_vec)
        if(~isempty(strmatch(upper(genome_version), upper(genome_versions_vec{i}{1}), 'exact'))) % enable also already to input the organism
            organism_str = genome_versions_vec{i}{1};
            return;
        end
        if(~isempty(strmatch(upper(genome_version), upper(genome_versions_vec{i}{2}), 'exact'))) % enable also already to input the organism
            organism_str = genome_versions_vec{i}{1};
            return;
        end
    end
    organism_str = 'HUMAN'; % default is human (if no version was found)
end




