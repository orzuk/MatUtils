if(read_vcf_flag)
    max_chr=23; min_chr=1;
else % run only once
    max_chr=1; min_chr=1;
end
sub_dir_str = dir_from_file_name(exome_struct.spectrum_data_file);
chr_file_str = suffix_from_file_name(remove_suffix_from_file_name(exome_struct.spectrum_data_file));

%                min_chr=22; max_chr=23; % TEMP FOR DEBUG! TAKE SHORT CHROMOSOME
for chr = min_chr:max_chr % Issue: not always good to split by chromosomes !! % take short chrom for debugging min_chr:max_chr % 1:23
    do_chr = chr
    if(read_vcf_flag) % One file per chromosome. (includes ALL populations)
        tmp_file_name = GetFileNames(fullfile(spectrum_data_dir, sub_dir_str, ...
            [exome_struct.prefix '*.chr' chr_num2str(chr) '.' chr_file_str '.vcf' ]));
        exome_struct.spectrum_data_file = fullfile(sub_dir_str, tmp_file_name{1});
        %                 spectrum_data_files{i} = fullfile(sub_dir_str, ...
        %                     ['ESP6500.chr' chr_num2str(chr) '.' chr_file_str '.vcf' ]); % '_' population{1} '.vcf']; %    .vcf'];
    else % One file. Already includes population string
        %                spectrum_data_files{i} = ['/Tennessen_Science_2012/all_chr_ESP6500.snps' '_' population{1} '.mat']; %    .vcf'];
        
        if(read_to_mat_flag)
            exome_struct.spectrum_data_file = fullfile(dir_from_file_name(exome_struct.spectrum_data_file), ...
                ['all_chr_' exome_struct.prefix '.' chr_file_str '.vcf']); %  '_' population{1} '.mat']; %    .vcf'];
        else % already have a .mat file
            exome_struct.spectrum_data_file = fullfile(dir_from_file_name(exome_struct.spectrum_data_file), ...
                ['all_chr_' exome_struct.prefix '.' chr_file_str '.mat']); %  '_' population{1} '.mat']; %    .vcf'];
        end
    end
    exome_struct.spectrum_data_files_str = [exome_struct.spectrum_data_files_str '''' ...
        remove_suffix_from_file_name(exome_struct.spectrum_data_file) '_' population{1} '.mat' ''','];
    
    if(parse_site_frequency_flag)
        job_str = ['[A] =' ... % , n_vec, count_vec, f_vec, allele_types] = ' ...
            'parse_site_frequency_data(''' fullfile(spectrum_data_dir, exome_struct.spectrum_data_file) ...
            ''', [], ' num2str(read_to_mat_flag) ', ' num2str(extract_fields_flag) ', ' ...
            num2str(compute_gene_matrices_flag) ');']; %, gene_list
        
        %                 [A n_vec count_vec f_vec allele_types] = ...
        %                     parse_site_frequency_data(fullfile(spectrum_data_dir, spectrum_data_files{i})); % , gene_list);
        
        if(unite_flag)
            %                    if(chr == min_chr)
            %                        A = load(fullfile(spectrum_data_dir, sub_dir_str, ...
            %                            ['all_chr_ESP6500.' chr_file_str '_'  population{1} '_up_to_chr' num2str(chr) '.mat'])); % load union
            %                    else % load current
            A = load([remove_suffix_from_file_name(fullfile(spectrum_data_dir, exome_struct.spectrum_data_file))  '_' population{1} '.mat']);
            A = my_rmfield(A, 'INFO_ARR');
            %                    end
            
            in_matlab_flag=1;
        else % don't unite chroms
            if(in_matlab_flag)
                eval(job_str);
            else
                SubmitMatlabJobToFarm(job_str, ...
                    fullfile('out', ['parse_ESP_chr' chr_num2str(chr) '.out']), queue_str, ...
                    [], [], mem_flag); % allow specifying memory allocation
            end
        end % if unite chroms
        if(read_vcf_flag && in_matlab_flag) % unite different files (even without 'unite chr')
            if(isfield(A, 'GENE'))
                A.GENE = vec2column(A.GENE);
            end
            if(unite_flag)
                if(chr == min_chr) % 1)
                    all_A = A;
                else
                    all_A = union_SFS_structs(all_A, A); 
                end
            end % unite chr
            close all;
        end % if read vcf
        if(unite_flag)  % Here unite - why only first population? (Europeans?)
            save(fullfile(spectrum_data_dir, sub_dir_str, ...
                ['all_chr_' exome_struct.prefix '.' chr_file_str '_'  population{1} '_up_to_chr' num2str(chr) '.mat']), '-struct', 'all_A'); % Save union
        end
    end % if parse SFS data
end % loop on chr.
