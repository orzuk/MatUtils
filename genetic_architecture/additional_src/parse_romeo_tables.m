% Parse supp info tables from Romeo et al. ICS 2009
% function parse_romeo_tables()
Assign24MAmmalsGlobalConstants;
supp_info_dir = '../../common_disease_model/data/rare/Romeo_2009_ICS_Tryglicerides/';
supp_tab_files = {'SuppTable1_ICS.txt', 'SuppTable2_ICS.txt', 'SuppTable3_ICS.txt'};


for i=1:length(supp_tab_files)
    R = loadcellfile(fullfile(supp_info_dir, supp_tab_files{i})); % load tab-delimited txt files 
    switch i
        case 1 % table with individual alleles
            header = empty_cell_to_empty_str(R(:,1));
            rs_inds = my_strmatch('rs', header); % deal with snp ids (currently we exclude them) 
            for j=1:length(rs_inds)
               header{rs_inds(j)} = -1;  
            end
            gene_start_inds = my_strmatch('ANGPTL', header); gene_stop_inds = [gene_start_inds(2:end)'-1  length(header)-1];
            AlleleStruct = [];  AlleleStruct.populations = {'European', 'African', 'Hispanic', 'Other'}; 
            AlleleStruct.num_individuals = [1045 1830 601 75]; % how many individuals in each population 
            AlleleStruct.strand = [POS_STRAND POS_STRAND REV_STRAND REV_STRAND]; % strand of all genes: ANGPTL3,4,5,6
            for j=1:length(gene_start_inds)
                AlleleStruct.gene_names{j} = header{gene_start_inds(j)}; 
                AlleleStruct.mutation{j} = R((gene_start_inds(j)+1):gene_stop_inds(j), 2); 
                AlleleStruct.mutation{j} = strrep_cell(AlleleStruct.mutation{j}, '>', '/'); % standard mutation symbol 
               AlleleStruct.type{j} = R((gene_start_inds(j)+1):gene_stop_inds(j), 3); 
               AlleleStruct.MAF{j} = cell2mat( R((gene_start_inds(j)+1):gene_stop_inds(j), 4:6) ); 
               AlleleStruct.counts{j} = bsxfun(@(x, y) round(x.*y./100), AlleleStruct.MAF{j}, 2.*AlleleStruct.num_individuals(1:3))
               AlleleStruct.positions{j} = cell2mat( header((gene_start_inds(j)+1):gene_stop_inds(j)) ); 
               AlleleStruct.chr(j) = str2nums(R{gene_start_inds(j),3});
               rare_inds = find(AlleleStruct.MAF{j}(:,1) < 0.01 * 100); % decide on frequency cutoff in Europeans
               AlleleStruct.Seqs{j} = ExtractSeqsByPositions(AlleleStruct.chr(j), ...
                   min(AlleleStruct.positions{j}(AlleleStruct.positions{j}>0)), max(AlleleStruct.positions{j})+1, ...
                    seqs_dir, 'hg18', 1); 
                
%                AlleleStruct.ExonSeqs = ExtractExons(
               missense_inds = strmatch('Nonsynonymous', AlleleStruct.type{j}); % Determine which non-synonymous are in fact frameshifts or stop codons 
               missense_ref_allele =  AlleleStruct.Seqs{j}.seqs{1}( AlleleStruct.positions{j}(missense_inds) - ...
                    AlleleStruct.Seqs{j}.pos_start_vec + 1 );
               tmp_c = strsplit_cell(AlleleStruct.mutation{j}, '/');
               for k=1:length(tmp_c)
                  AlleleStruct.RefAllele{j}{k} = tmp_c{k}{1};
                  AlleleStruct.DerAllele{j}{k} = tmp_c{k}{2};                  
               end
               for k=1:length(missense_inds) % get codons (need to get right both strand and 3-split !!!
                   split = 1;
                   AlleleStruct.RefCodon{j}{k} = '';
                   switch AlleleStruct.strand(j) % get strand to get the codon
                       case POS_STRAND
                           AlleleStruct.DerCodon{j}{k} = AlleleStruct.RefCodon{j}{k};
                           AlleleStruct.DerCodon{j}{k}(split) = 'A'; % AlleleStruct.DerAllele{j}{missense_inds(k)};
                       case REV_STRAND
                           AlleleStruct.RefCodon{j}{k} = seqrcomplement(AlleleStruct.RefCodon{j}{k}); % use matlab toolbox
                           AlleleStruct.DerCodon{j}{k} = 'A'; % TEMP!
                   end
               end
               
               frameshift_inds = strmatch('del', AlleleStruct.DerAllele{j}); 
               stop_inds = []; %strmatch('',  AlleleStruct.DerCodon{j}); % find stop-codons 
               
               missense_rare_inds = setdiff( intersect(rare_inds, missense_inds), frameshift_inds ); 
               nonsense_rare_inds = intersect(rare_inds, frameshift_inds);
               AlleleStruct.missense_cum_freq(j) = sum(AlleleStruct.MAF{j}(missense_rare_inds));
               AlleleStruct.missense_cum_counts(j) = sum(AlleleStruct.counts{j}(missense_rare_inds)); 
               AlleleStruct.nonsense_cum_freq(j) = sum(AlleleStruct.MAF{j}(nonsense_rare_inds));
               AlleleStruct.nonsense_cum_counts(j) = sum(AlleleStruct.counts{j}(nonsense_rare_inds));
               
            end % loop on gene 
            save(file_name_to_mat(fullfile(supp_info_dir, supp_tab_files{i})), 'AlleleStruct');             
            
        case 2 % summary statistics for all alleles
            
        case 3 % summary statistics for phenotype
            
            
    end
    
    
end
    