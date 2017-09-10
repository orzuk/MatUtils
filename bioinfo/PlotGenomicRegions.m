% Plot the conservation of each type of genomic region 
%
% Input: 
% regions_file - file with regions information
% k - kmer size
% iters - # of iterations for each region type
% alpha - top conservation level (e.g. 5%)
% genome_version - what genome to work on
% fig_outfile - where to save the figure 
%
function Dummy = PlotGenomicRegions(regions_file, k, iters, alpha, genome_version, fig_outfile, varargin)

Assign24MammalsGlobalConstants();

R = load(regions_file);  % load regions 
upstream = 2500; downstream = 0;
Promoters = ExtractPromoters('../data', genome_version, upstream, downstream, 1, [], 1);

R.chr_vec = [vec2row(R.chr_vec) vec2row(Promoters.chr_vec)]';
R.pos_start_vec = [vec2row(R.pos_start_vec) vec2row(Promoters.pos_start_vec)]';
R.pos_end_vec = [vec2row(R.pos_end_vec) vec2row(Promoters.pos_end_vec)]';
R.strand_vec = [vec2row(R.strand_vec) vec2row(Promoters.strand)]';
R.gene_names_vec = [vec2row(R.gene_names_vec) vec2row(Promoters.gene_names)]';
R.region_type_vec = [vec2row(R.region_type_vec) repmat(PROMOTER, 1, length(Promoters.chr_vec))]';

unique_region_types = unique(R.region_type_vec); 
num_region_types = length(unique_region_types);

R.region_lengths = R.pos_end_vec - R.pos_start_vec + 1; % get region lengths

region_types_lengths = zeros(num_region_types,1); 
for i=1:num_region_types % loop over region types
    cur_reg_inds = find(R.region_type_vec == unique_region_types(i));
    region_types_lengths(i) = sum(R.region_lengths(cur_reg_inds))
end

Dummy = 0;

for i=1:num_region_types % loop over region types
    cur_reg_inds = find(R.region_type_vec == unique_region_types(i));
    rand_regions = GetRandomGenomicRegions(iters, k, genome_version, ...
        R.chr_vec(cur_reg_inds), R.pos_start_vec(cur_reg_inds), R.pos_end_vec(cur_reg_inds)); % get some random regions from the genome (NOT input regions)

    pwms_flag = 0; loglike_flag = 1; positions_flag = 0; tree_flag = 0; branch_length_flag = 0;
    PI = ExtractPWMsByPositions(rand_regions.chr_vec, rand_regions.pos_start_vec, rand_regions.pos_end_vec, ...
        [], genome_version, ...
        pwms_flag, loglike_flag, positions_flag, tree_flag, branch_length_flag);
    kmers_pi{i} = sum_cell(PI.loglike); 
end


%%%%%%%%% Plot cumulative curve of conservation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; hold on;
for i=1:num_region_types % loop over region types
    n = length(kmers_pi{i});
    plot(sort(kmers_pi{i}), [1:n]./n, color_vec(i)); 
end
h_legend = legend(genome_types, 4);
xlabel('\Pi LODs'); ylabel('cum. freq.');
if(exist('fig_outfile', 'var'))
    saveas(gcf, [fig_outfile '.regions_conservation_dist.fig']);
    saveas(gcf, [fig_outfile '.regions_conservation_dist.jpg']);
end
%%%%%%%%%%%%%%%% Now without the non-aligned parts 
figure; hold on;
for i=1:num_region_types % loop over region types
    pos_kmers = kmers_pi{i}(kmers_pi{i} > 0);
    n = length(pos_kmers);
    plot(sort(pos_kmers), [1:n]./n, color_vec(i)); 
end
h_legend = legend(genome_types, 4);
xlabel('\Pi LODs'); ylabel('cum. freq.');
if(exist('fig_outfile', 'var'))
    saveas(gcf, [fig_outfile '.regions_conservation_dist_pos.fig']);
    saveas(gcf, [fig_outfile '.regions_conservation_dist_pos.jpg']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%% Pie chart of Genome composition: All genome vs. conserved parts %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;  subplot(1,2,1);  % numeric values to different genome regions
pie(region_types_lengths); title('genome'); h_legend = legend(genome_types, 3);  set(h_legend,'FontSize',7);

save('temp_kmers_pi.mat', 'PI', 'kmers_pi');

region_types_frac_conserved = zeros(num_region_types, 1); 
q = quantile(cell2vec(kmers_pi), 1-alpha);
for i=1:num_region_types % loop over region types
   region_types_frac_conserved(i) = sum(kmers_pi{i} >= q) ./ iters; 
end
subplot(1,2,2);  % numeric values to different genome regions - conserved !!!
pie(region_types_lengths .* region_types_frac_conserved); title('conserved'); h_legend = legend(genome_types, 3);  set(h_legend,'FontSize',7);
if(exist('fig_outfile', 'var'))
    saveas(gcf, [fig_outfile '.regions_conservation_pie.fig']);
    saveas(gcf, [fig_outfile '.regions_conservation_pie.jpg']);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% check overlap ... 
% Extract kmers 

Dummy = 0;




