% Load all the probesets of a chip
function probesets = get_probesets(chip)

%load(fullfile('..','database', [chip '_probes.mat']));
genome_assembly = get_genome_assembly();

load(['..\database\' chip '_annot_data_' genome_assembly '.mat'], 'snp_ids');
probesets = snp_ids;