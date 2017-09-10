% The function loads from the database the ending of the p-locations
% of all the chromosomes
function end_p_location = load_end_p_loc()

%% load('chr_data.mat', 'end_p_location'); % Old: Version without the database

genome_assembly = get_genome_assembly();
load(['../Database/chr_data_' genome_assembly '.mat'], 'end_p_location');  % New: Version with the database
