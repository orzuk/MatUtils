% Written by Or Zuk 7/2007
%
% This script/function Performs all the updates needed for the Hapmap data
% once a new version is available. It is assumed that all the .txt files
% were already downloaded and put in the right places.
% Note: It is important when calling this function that the hapmap_version 
% will agree with the hapmap_data_http. To ensure this, we set the
% convention that the hapmap_version must be a substring of the file names 
% that we download. 
%
% The Input: 
% hapmap_version - current version of Hapmap 
% hapmap_dir - directory for Hapmap files
% chip_types - which types of chips to support
%
function Dummy  = UpdateDatabaseHapmap(hapmap_dir, chip_types)
AssignAllGlobalConstants;

save_flag = 1; % We need to save new files
chroms = 1:24; % include X and Y chromosomes

% get correct hapmap version 
hapmap_version = genome_assembly_to_hapmap_version()

% get correct url
hapmap_data_http = hapmap_version_to_http(hapmap_version);




file_names = get_hapmap_file_names_from_web(hapmap_version);

mat_file_names = file_names;

for i=1:length(mat_file_names)
    txt_idx = strfind(file_names{i}, 'txt');
    mat_file_names{i} = [file_names{i}(1:txt_idx-1) '.mat'];
end

% s = urlread(hapmap_data_http); % Read data from project
% start_inds = strfind(s,'genotypes_chr'); end_inds = strfind(s,'txt.gz')+5;
% if(size(start_inds) ~= size(end_inds))
%     sprintf('Error In Page - Go and Download it Manually ...')
%     Dummy = [];
%     return;
% end
% 
% file_names = {};
% for i=1:length(start_inds)
%     file_names{i} = s(start_inds(i):end_inds(i));
% end
% file_names = unique(file_names);
% 
% if(isempty(strfind(file_names{1},hapmap_version)))
%     sprintf('Error In Hapmap Version - Go and Download it Manually ...')
%     Dummy = [];
%     return;
% end

for hapmap_population = [CEU JPT_CHB YRI]
    download_pop = hapmap_population

    % First load the data and extract it to the correct places
    download_ctr=0;
    for i=1:length(file_names)
        if(~isempty(findstr(pop_str_vec_plus{hapmap_population}, file_names{i})))
            download_ctr=download_ctr+1
            
            if( (~exist( fullfile(hapmap_dir, pop_str_vec{hapmap_population}, hapmap_version, file_names{i}(1:end-3)) )) && ...
                (~exist( fullfile(hapmap_dir, pop_str_vec{hapmap_population}, hapmap_version, mat_file_names{i}) ))   )               % don't waste time downloading again and again
                gunzip([hapmap_data_http, file_names{i}],fullfile(hapmap_dir, pop_str_vec{hapmap_population}, hapmap_version));
            end
        end
    end
end

for hapmap_population = [CEU JPT_CHB YRI] 
    pop = hapmap_population
    ctr=0;
    for i=1:length(file_names)
        if(~isempty(findstr(pop_str_vec_plus{hapmap_population}, file_names{i})))
            ctr=ctr+1;
            pop_file_name.txt{ctr} = file_names{i}(1:end-3); % remove the '.gz' (gunzip) at the end 
            pop_file_name.mat{ctr} = [file_names{i}(1:end-6) 'mat']; % remove the '.gz' (gunzip) at the end 
        end
    end    
    
    [SnpsData, SnpsBases, SnpsFreqs, SnpsHetros, SnpsBadCalls, SnpsChromLocs SnpsNames] = ...
        HapMapGenotypeTxtToMatFiles(hapmap_population, hapmap_version, hapmap_dir, chroms, pop_file_name.txt, save_flag); % Convert all .txt files to .mat files

    for i=1:length(chip_types)  % Intersect with all current Chips
        chip = chip_types{i}
        Dummy = HapMapExtractOnlyChipGenotypes(chroms, chip_types{i}, hapmap_version, hapmap_dir, hapmap_population, pop_file_name.mat);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This part is used for updating the LD structure
use_hetro_flag = 0; % We currently don't know how to use them
for hapmap_population = [CEU] %  JPT_CHB YRI] % NOTE: Currently use only Europeans. Should add also others !!! 
    for i=1:length(chip_types)  % Intersect with all current Chips
        chip = chip_types{i}
        LD_dummy = ComputeSNPsAdjacentLD(chroms, hapmap_dir, chip, hapmap_population, hapmap_version, genome_assembly, use_hetro_flag); % this function saves the LD structures
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

