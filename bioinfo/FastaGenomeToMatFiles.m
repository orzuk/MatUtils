% Convert whole genome chromsomes in fasta format into binary .mat files. 
% Current genome version 'supported' are Hg17 and Hg18. The positions are just
% consecutive and it should work fine (but we need to check this ....)
% 
% Input: 
% genome_dir - directory with fasta genome files
% genome_version - which genome version do we use
% output_dir - where to write the output .mat files to
% organism_str - which organism is it
function Dummy = FastaGenomeToMatFiles(genome_dir, genome_version, output_dir, organism_str)

Assign24MammalsGlobalConstants;

if(~exist('organism_str', 'var')) % default organism should be human
    organism_str = 'HUMAN';
end

num_chroms = organism_chr; 
chroms = [1:num_chroms+2];  % take also mitochondria .. change back to one !!!

% hg_dir = ['/seq/genome/human/human_hg' num2str(genome_ver)];
% data_dir = ['../data/hg' num2str(genome_ver)];
if(strcmp(upper(organism_str), 'HUMAN'))
    GenomeVersion = ['hg' num2str(genome_version)];
end
if(strcmp(upper(organism_str), 'MOUSE'))
    GenomeVersion = ['mm' num2str(genome_version)];
end

for i=1:length(chroms)
    chr = chroms(i)
    if(exist(fullfile(genome_dir, chr_num2str(chr, organism_str)), 'dir'))
        [Header Sequence] = fastaread(fullfile(genome_dir, chr_num2str(chr, organism_str), ['chr' chr_num2str(chr, organism_str) '.fa']));
    else
        [Header Sequence] = fastaread(fullfile(genome_dir, ['chr' chr_num2str(chr, organism_str) '.fa']));
    end
    %    load(fullfile(data_dir, ['chr' chr_num2str(chr) '.mat'] ));
    save(fullfile(output_dir, ['chr' chr_num2str(chr, organism_str) '_case_sensitive.mat'] ), 'Header', 'GenomeVersion', 'Sequence');
    Sequence = upper(Sequence); % get rid of size ...
    save(fullfile(output_dir, ['chr' chr_num2str(chr, organism_str) '.mat'] ), 'Header', 'GenomeVersion', 'Sequence');
end

Dummy = 0;



