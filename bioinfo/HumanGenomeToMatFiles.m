% Convert all the human genome into binary .mat files. 
% Current genome version 'supported' are Hg17 and Hg18. The positions are just
% consecutive and it should work fine (but we need to check this ....)
function Dummy = HumanGenomeToMatFiles(genome_ver)

% Assign24MammalsGlobalConstants;
hg_dir = ['/seq/genome/human/human_hg' num2str(genome_ver)];
data_dir = ['../data/hg' num2str(genome_ver)];
chroms = 1:25;  % take also mitochondria .. change back to one !!!
GenomeVersion = ['hg' num2str(genome_ver)];

for i=1:length(chroms)
    chr = chroms(i)
    read_input_file = fullfile(hg_dir, chr_num2str(chr), ['chr' chr_num2str(chr) '.fa'])
    [Header Sequence] = fastaread(fullfile(hg_dir, chr_num2str(chr), ['chr' chr_num2str(chr) '.fa']));
    %    load(fullfile(data_dir, ['chr' chr_num2str(chr) '.mat'] ));
    save(fullfile(data_dir, ['chr' chr_num2str(chr) '_case_sensitive.mat'] ), 'Header', 'GenomeVersion', 'Sequence');
    Sequence = upper(Sequence); % get rid of size ...
    save(fullfile(data_dir, ['chr' chr_num2str(chr) '.mat'] ), 'Header', 'GenomeVersion', 'Sequence');
end

Dummy = 0;



