% Compute some statistics on the genomic regions
% 
% Input: 
% genomic_regions - a struct containing (chr start stop) or a file name containing regions
%
% Output:   
% regions_total_lens - total length of all regions
% regions_nums - number of regions
% regions_mean_lens - mean of regions lens
% regions_std_lens - std of regions lens
% regions_names - names of regions
%
function [regions_total_lens regions_nums regions_mean_lens regions_std_lens regions_names] = ...
    GenomicRegionsStatistics(genomic_regions)

doing_types = 0;
if(doing_types)    % old stuff: get the total length of each genomic type
    regions_total_lens = zeros(1,5);
    for i=1:5
        regions_total_lens(i) = sum( genomic_regions.pos_end_vec(genomic_regions.region_type_vec == i) - ...
            genomic_regions.pos_start_vec(genomic_regions.region_type_vec == i) );
    end
else
    if(ischar(genomic_regions)) % this means that we're given an input file/dir
        params_dir = genomic_regions;
        chip_seq_params_files = setdiff(GetFileNames(fullfile(params_dir, '*chip*.txt')), ...
        GetFileNames(fullfile(params_dir, '#*chip*'))); % avoid some left-over files from unix 
        num_regions_sets = length(chip_seq_params_files);
        ctr=1;
        for i=1:num_regions_sets
            i_is = i
            motif_params = ReadMotifFindParametersFile(fullfile(params_dir, chip_seq_params_files{i})); % read parameters
            clear tf_names;
            load(fullfile(motif_params.regions_dir, motif_params.regions_file));
            if(iscell(chr_vec))
                for j=1:length(chr_vec)
                    genomic_regions.chr_vec{ctr} = chr_vec{j};
                    genomic_regions.pos_start_vec{ctr} = pos_start_vec{j};
                    genomic_regions.pos_end_vec{ctr} = pos_end_vec{j};
                    genomic_regions.name{ctr} = tf_names{ctr}; 
                    ctr=ctr+1;
                end
            else
                genomic_regions.chr_vec{ctr} = chr_vec;
                genomic_regions.pos_start_vec{ctr} = pos_start_vec;
                genomic_regions.pos_end_vec{ctr} = pos_end_vec;
                if(exist('tf_names', 'var'))
                    if(iscell(tf_names))
                        genomic_regions.name{ctr} = tf_names{1};
                    else
                        genomic_regions.name{ctr} = tf_names;
                    end
                else
                    genomic_regions.name{ctr} = motif_params.regions_file;
                end
                ctr=ctr+1;
            end
        end
    end
    if(iscell(genomic_regions.pos_start_vec))           
        num_regions_sets = length(genomic_regions.pos_start_vec);
    else
        num_regions_sets = 1;
    end
    regions_total_lens = zeros(1, num_regions_sets);
    regions_nums = zeros(1, num_regions_sets);
    regions_mean_lens = zeros(1, num_regions_sets);
    regions_std_lens = zeros(1, num_regions_sets);
    regions_names = cell(1, num_regions_sets); 
    for i=1:num_regions_sets % loop over different sets of regions
        if(iscell(genomic_regions.pos_start_vec))
        regions_lens_vec = genomic_regions.pos_end_vec{i} - genomic_regions.pos_start_vec{i};  % get the vector of regions lens
        else
            regions_lens_vec = genomic_regions.pos_end_vec - genomic_regions.pos_start_vec;  % get the vector of regions lens
        end
        regions_total_lens(i) = sum(regions_lens_vec);
        regions_nums(i) = length(regions_lens_vec);
        regions_mean_lens(i) = mean(regions_lens_vec);
        regions_std_lens(i) = std(regions_lens_vec);
        if(isfield(genomic_regions, 'name'))
            if(iscell(genomic_regions.name))
                regions_names{i} = genomic_regions.name{i};
            else
                regions_names{i} = genomic_regions.name{i};
            end
        else
            if(isfield(genomic_regions, 'tf_names'))
                if(iscell(genomic_regions.tf_names))
                    regions_names{i} = genomic_regions.tf_names{i};
                else
                    regions_names{i} = genomic_regions.tf_names;
                end
            end
        end
        figure; hist(regions_lens_vec, 100); title(['regions for ' regions_names{i}]);
    end
end