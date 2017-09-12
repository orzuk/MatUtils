%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Choose operations to run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AssignGeneralConstants;
run_examples = 1;
run_examples_str = 'david'; % David % 'hartl'; % 'new_runs' % which simulations to run

concatenate_chr = 0;
run_preprocess = 0;
preprocess_schizophrenia = 0;
prepare_LD_matrices = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ttt = cputime;

eliana_SNPs_data_dir = '/broad/landerlab2/eliana/Pipeline/SNPTEST_OUTPUT'
SNPs_data_dir = '../../common_disease_model/data/fdr';
SNPs_figs_dir = fullfile(SNPs_data_dir, 'figs');
LD_file = fullfile(SNPs_data_dir, 'PLINK', 'CEU_chr1.ld_preprocessed.mat');
chr1_LD_file = fullfile(SNPs_data_dir, 'PLINK', 'CEU_chr1.ld_preprocessed.mat');


genome_LD_file = fullfile(SNPs_data_dir, 'PLINK', 'CEU_genomic.ld_preprocessed.mat'); % combined all chromosomes
causal_SNPs_file = fullfile(SNPs_data_dir, 'causal_loci_chr1.txt');

%GWAS_file = fullfile(SNPs_data_dir, 'chr1_20_causal_snps_bin_40_50.txt'); % always include a suffix, otherwise things get screwed up
GWAS_file = fullfile(SNPs_data_dir, 'chr1_20_causal_snps_bin_40_50_only_causal_no_LD_with_FDR.txt'); % always include a suffix, otherwise things get screwed up
schizophrenia_GWAS_file = fullfile('schizophrenia', 'SCZ_data_wfdr'); % no dir, no suffix
bipolar_GWAS_file = fullfile('bipolar', 'BIP_data_wfdr'); % bipolar
height_GWAS_file = fullfile('height', 'GIANT_data_wfdr'); % height
% prune_SNPs_by_LD_script([],[],bipolar_GWAS_file,[],0);

GWAS_file_str = {'chr1_20_causal_snps_bin_40_50_only_causal_no_LD_with_FDR', ...
    'sample_eliana/twenty_uniform_new.out_columns', 'sample_eliana/twenty_lowfreq_new.out_columns', ...
    schizophrenia_GWAS_file, height_GWAS_file, bipolar_GWAS_file, ...
    'sample_eliana/all_lower_0.1_upper_0.2_h2_0.8_COLUMNS_10000_WHEADER', ...
    'sample_eliana/all_lower_0.1_upper_0.9_h2_0.8_COLUMNS_10000_WHEADER', ...
    'sample_eliana/all_lower_0.4_upper_0.5_h2_0.8_COLUMNS_10000_WHEADER', ...
    'sample_eliana/all_lower_0.1_upper_0.2_h2_0.8_ALLPHENO_COLUMNS_10000_WHEADER', ...
    'sample_eliana/all_lower_0.1_upper_0.9_h2_0.8_ALLPHENO_COLUMNS_10000_WHEADER', ...
    'sample_eliana/all_lower_0.4_upper_0.5_h2_0.8_ALLPHENO_COLUMNS_10000_WHEADER', ... % start for new
    'sample_eliana/all_lower_0.01_upper_0.02_h2_0.8_columns_10000_wheader', ...
    'sample_eliana/all_M1_0.01_M2_0.4_h2_0.8_columns_10000_wheader', ...
    'sample_eliana/all_M1_0.2_M2_0.4_h2_0.8_columns_10000_wheader'}; %  height_GWAS_file}; % possible file names for GWAS results
trait_type_str = {'quantitative', 'quantitative', 'quantitative', ...
    'disease', 'quantitative', 'disease', ...
    'quantitative', 'quantitative', 'quantitative', ...
    'quantitative', 'quantitative', 'quantitative', ...
    'quantitative', 'quantitative', 'quantitative'};

if(run_examples)
    perform_deconvolution = 0; % 1; % flag saying if to perform deconvolution
    %    prune_SNPs_by_LD_script([],LD_file,GWAS_file,causal_SNPs_file,0); % a file with 20 loci in the highest bin (40-50%)
    
    rand_GWAS_file = fullfile(SNPs_data_dir, 'randresults_foror');
    %prune_SNPs_by_LD_script([],LD_file,rand_GWAS_file,[],0);
    
    % prune_SNPs_by_LD_script([],[],schizophrenia_GWAS_file,[],0);
    
    all_GWAS_files = GetFileNames(fullfile(SNPs_data_dir, ...
        'sample_eliana', run_examples_str, 'all_*beta*columns_*.txt'), 1); % New! choose files (Hartl/Peaks)
    
    if(prepare_LD_matrices)
        prune_SNPs_by_LD_script([],genome_LD_file,GWAS_file,causal_SNPs_file, ...
            'prepare_LD_matrices', 0, trait_type_str{i});        % run all. Use genomic LD file. Don't save as .txt
    end % run over all chromosomes and prepare LD matrics
    for i=1:length(all_GWAS_files)
        %        13:15 % 7:length(GWAS_file_str) % 2: % loop on different simulations
        run_gwas_results_file = i
        GWAS_file = all_GWAS_files{i}; % GWAS_file = fullfile(SNPs_data_dir, [GWAS_file_str{i} '.txt']);
        tmp_causal_ind = strfind(lower(GWAS_file), 'columns');
        tmp_end_ind = find(GWAS_file == '_');
        tmp_end_ind = tmp_end_ind(find(tmp_end_ind > tmp_causal_ind, 2));
        switch run_examples_str
            case 'david'
                causal_SNPs_file = [];
                SNPs_file = fullfile(SNPs_data_dir, 'sample_eliana', 'all_snp_freqs.txt');
            otherwise
                causal_SNPs_file = strrep(GWAS_file, GWAS_file(tmp_causal_ind:tmp_end_ind(2)), 'diseaseloci_');
                SNPs_file = []; % no need (we've got the allele frequencies) 
        end
        trait_type_str{i} = 'quantitative';
        %         switch i
        %             case {1,2,3} % simulated data
        %                 causal_SNPs_file = fullfile(SNPs_data_dir, [GWAS_file_str{i} '.causal.txt']);
        %             case {7,8,9,10,11,12,13,14,15} % new simulatied data
        %                 tmp_causal_str = strsplit('.', GWAS_file_str{i});
        %                 tmp_causal_str = [cell2vec(tmp_causal_str(1:end-1), '.') '.' tmp_causal_str{end}(1)]
        %                 causal_SNPs_file = fullfile(SNPs_data_dir, [tmp_causal_str '_diseaseloci.txt']);
        %             otherwise % real data
        %                 causal_SNPs_file = [];
        %         end
        switch machine
            case PC
                prune_SNPs_by_LD_script(SNPs_file,genome_LD_file,GWAS_file,causal_SNPs_file, ...
                    0,0, trait_type_str{i}, perform_deconvolution);        % run all. Use genomic LD file. Don't save as .txt
            case UNIX
                eval_str = ['prune_SNPs_by_LD_script([],''' genome_LD_file ''',''' GWAS_file ''',''' causal_SNPs_file ''', ' ...
                    '0,0, ''quantitative'' , ' num2str(perform_deconvolution) ');'];         % run all. Use genomic LD file. Don't save as .txt
                SubmitMatlabJobToFarm(eval_str, ...
                    ['out/run_LD_prunning_' remove_dir_from_file_name(GWAS_file) '.out'], 'hour');
        end % switch machine
        
        
        close all; % not to put to many figures
    end % loop on GWAS files 
end % run examples

if(concatenate_chr) % Cimbine pairwise LD data from different chromosomes 
    R_LD = []; R_LD.inds = [];   R_LD.rsid = [];  R_LD.chr = [];  R_LD.pos = [];  R_LD.LD_inds = [];  R_LD.R2 = [];
    counter = 0;
    for chr=1:22 % unite LD blocks
        concatenate_chr = chr
        chr_LD_file = fullfile(SNPs_data_dir, 'PLINK', ['CEU_chr' num2str(chr) '.ld_preprocessed.mat']);
        CHR = load(chr_LD_file);
        R_LD.inds = [R_LD.inds' CHR.R_LD.inds'+counter]';        
        R_LD.rsid = [R_LD.rsid' CHR.R_LD.rsid']';
        R_LD.chr = [R_LD.chr' CHR.R_LD.chr']';
        R_LD.pos = [R_LD.pos' CHR.R_LD.pos']';
        R_LD.LD_inds = [R_LD.LD_inds' sum_cell(CHR.R_LD.LD_inds', counter)]';
        R_LD.R2 = [R_LD.R2' CHR.R_LD.R2']';
        counter = counter + length(CHR.R_LD.inds); % update counter: cumulative # of SNPs
    end % loop on chromosomes
    save(fullfile(SNPs_data_dir, 'PLINK', ['CEU_genomic.ld_preprocessed.mat']), 'R_LD');
    WriteDataFile(R_LD, fullfile(SNPs_data_dir, 'PLINK', ['CEU_genomic.ld_preprocessed.txt']));
end % if concatenate



if(run_preprocess) % Preprocess pairwise LD data in each chromosome
    switch machine % concatenate LD information for all SNPs
        case {UNIX, PC} % why only in UNIX ???
            for chr = 1:22 % Run all autosomal chromosomes preprocessing
                run_LD_preprocess_chr = chr
                chr_LD_file = fullfile(SNPs_data_dir, 'PLINK', ['CEU_chr' num2str(chr) '.ld.txt']);
                chr_job_str = ['prune_SNPs_by_LD_script([],''' chr_LD_file ''',[],[],1)'];
                switch machine
                    case UNIX
                        SubmitMatlabJobToFarm(chr_job_str, ['out/LD_preprocess_chr' num2str(chr) '.out'], 'priority');
                    case PC
                        eval(chr_job_str);
                end
            end
    end
end % run preprocess

run_LD_time = cputime - ttt


if(preprocess_schizophrenia) % remove bad SNPS
    S = ReadDataFile('../../common_disease_model/data/fdr/schizophrenia/daner_SCZ17f_bad_snps_rsid.txt'); % Special preprocessing for scizophrenia
    load('../../common_disease_model/data/fdr/schizophrenia/SCZ_data_wfdr.mat'); % load GWAS_RES. Special preprocessing for scizophrenia
    [bad_snps bad_I bad_J] = intersect(GWAS_RES.rsid, S.SNP);
    good_I = setdiff(1:length(GWAS_RES.rsid), bad_I);
    GWAS_RES = struct_by_inds(GWAS_RES, good_I);
    save('../../common_disease_model/data/fdr/schizophrenia/SCZ_data_wfdr.mat', 'GWAS_RES');
end




