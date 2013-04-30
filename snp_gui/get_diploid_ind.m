%function diploid_ind = get_diploid_ind(sample_names)
function diploid_ind = get_diploid_ind(sample_names)

num_samples = length(sample_names);

diploid_ind = [];
for i = 1:num_samples
    sample_name = char(sample_names{i});
    if(strcmp(sample_name(end-1:end), '_n'))
        diploid_ind = [diploid_ind i];
    end
end