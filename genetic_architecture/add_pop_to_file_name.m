% Add population to SFS file name
function pop_file_name = add_pop_to_file_name(site_frequency_file_name, population)

pop_file_name = fullfile(dir_from_file_name(site_frequency_file_name), strdiff(population, '_'), ...
    [remove_suffix_from_file_name(remove_dir_from_file_name(site_frequency_file_name)) population '.mat']);

