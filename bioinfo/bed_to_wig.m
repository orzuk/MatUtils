% Convert .bed format (each element is transcript) to .wig (each element an exon)
function [chr_vec pos_start_vec pos_end_vec data] = bed_to_wig(bed_data)

total_exons = sum(bed_data.num_exons); % spread all the relative bed file to chr start end
chr_vec = zeros(total_exons,1);
pos_start_vec = zeros(total_exons,1);
pos_end_vec = zeros(total_exons,1);
ctr = 1;
for i=1:length(bed_data.num_exons)
    chr_vec(ctr:ctr+bed_data.num_exons(i)-1) = bed_data.chr_vec(i);
    for j=1:bed_data.num_exons(i) % so far ignore strand
        if(bed_data.strand{i} == '+')
            pos_start_vec(ctr+j-1) = ...
                bed_data.pos_start_vec(i) + bed_data.exon_starts{i}(j);
            pos_end_vec(ctr+j-1) = ...
                pos_start_vec(ctr+j-1) + bed_data.exon_lengths{i}(j);
        else % nagative strnad
            %       exons.pos_start_vec(ctr+j-1) = LINCS.pos_start_vec(i) + bed_data.exon_starts{i}(j);
            %       exons.pos_end_vec(ctr+j-1) = exons.pos_start_vec(ctr+j-1) + bed_data.exon_lengths{i}(j);
        end
    end
    if(bed_data.strand{i} == '+')
        ctr = ctr + bed_data.num_exons(i);
    end
end
chr_vec = chr_vec(1:ctr-1);
pos_start_vec = pos_start_vec(1:ctr-1);
pos_end_vec = pos_end_vec(1:ctr-1);


data = bed_data.gene_names; 
