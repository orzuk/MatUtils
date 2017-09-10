%function [sample_pairs_cell, pairs_samples_ind non_paired_ind_n non_paired_ind_d] = samples_into_pairs_old(samples);
% rem in first column, diag in second
function [sample_pairs_cell, pairs_samples_ind non_paired_ind_n non_paired_ind_d] = samples_into_pairs_old(samples);

if(size(samples, 1) ~= 1)
    samples = samples';
end

num_samples = length(samples);
diag_samples = cell(1,0);
diag_samples_ind = 0;
rem_samples = cell(1,0);
rem_samples_ind = 0;
non_paired_ind_n = [];
non_paired_ind_d = [];

%arrange in two groups
for i = 1:num_samples
    sample = char(samples(1,i));
    underscore_ind = strfind(sample, '_');
    t=strcmp(sample(underscore_ind(end)+1), 'n');
    if(t==1)
        rem_samples_ind = rem_samples_ind+1;
        rem_samples{1, rem_samples_ind} = sample(1:underscore_ind(end)-1);
        non_paired_ind_n = [non_paired_ind_n i];
    end
    t=strcmp(sample(underscore_ind(end)+1), 'd');
    if(t==1)
        diag_samples_ind = diag_samples_ind+1;
        diag_samples{1, diag_samples_ind} = sample(1:underscore_ind(end)-1);
        non_paired_ind_d = [non_paired_ind_d i];
    end
end

samples_with_pairs = intersect(rem_samples, diag_samples);
num_pairs = length(samples_with_pairs);

sample_pairs_cell = [];
if(num_pairs>0)
    sample_pairs_cell = cell(num_pairs, 2);
    sample_pairs_cell(:,1) = add_suff_to_cell_str(samples_with_pairs', '_n');
    sample_pairs_cell(:,2) = add_suff_to_cell_str(samples_with_pairs', '_d');
end

pairs_samples_ind = zeros(size(sample_pairs_cell));
for i = 1:num_pairs
    t = strmatch(sample_pairs_cell{i,1}, samples);
    pairs_samples_ind(i,1) = t;
    t = strmatch(sample_pairs_cell{i,2}, samples);
    pairs_samples_ind(i,2) = t;
end
if(num_pairs)
    non_paired_ind_d = setdiff(non_paired_ind_d, pairs_samples_ind(:,2)');
    non_paired_ind_n = setdiff(non_paired_ind_n, pairs_samples_ind(:,1)');
end
