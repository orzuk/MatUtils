function err=make_ucsc_files(ucsc_mounted_folder, segments_del, mat_genes_del, segments_samples_del, segments_samples_arm_del, ...
    segments_amp, mat_genes_amp, segments_samples_amp, segments_samples_arm_amp, ...
    genes_db_struct, chip_annot, sample_names, normal_disease_suffix)

err='';
[sample_names idx]=sort(sample_names);
segments_samples_arm_del=segments_samples_arm_del(:,idx);
segments_samples_arm_amp=segments_samples_arm_amp(:,idx);

genes_unique=union(mat_genes_del(:,1),mat_genes_amp(:,1));
for i=1:length(genes_unique)
    idx=find(strcmp(genes_db_struct.gene_symbols,genes_unique{i}));
    region_start=genes_db_struct.loc_start(idx(1))-1000000;
    if region_start<0
        region_start=0;
    end
    region_end=genes_db_struct.loc_end(idx(1))+1000000;
    chr=genes_db_struct.chr(idx(1));
    
    %find only segments overlapping in this region
    idx=find( segments_del(:,4) == chr & ...
               ( ( segments_del(:,1)<=region_start & segments_del(:,2)>=region_end ) | ...
                (segments_del(:,1)>= region_start & segments_del(:,1)< region_end  ) | ...
                (segments_del(:,2)> region_start & segments_del(:,2)<= region_end  ) ) );
    segments_del_in_region=segments_del(idx,:);
    segments_samples_del_in_region=segments_samples_del(idx);
    segments_samples_arm_del_in_region=segments_samples_arm_del(idx,:);
    
    idx=find(   segments_amp(:,4) == chr & ...
                ( ( segments_amp(:,1)<=region_start & segments_amp(:,2)>=region_end ) | ...
                (segments_amp(:,1)>= region_start & segments_amp(:,1)< region_end  ) | ...
                (segments_amp(:,2)> region_start & segments_amp(:,2)<= region_end  ) ) );
    segments_amp_in_region=segments_amp(idx,:);
    segments_samples_amp_in_region=segments_samples_amp(idx);
    segments_samples_arm_amp_in_region=segments_samples_arm_amp(idx,:);
            
    str1=['browser position chr' num2str(chr) ':' num2str(region_start) '-' num2str(region_end)];
    for j=1:length(sample_names)
        str1=[str1 '\n' 'track name=' sample_names{j} ' itemRgb="On" visibility="dense"'];
        
        if ~isempty(find(strcmp([segments_samples_arm_del_in_region(:,j); segments_samples_arm_amp_in_region(:,j)],'DEL')))
            str1=[str1 '\n' 'chr' num2str(chr) ' ' num2str(region_start) ' ' ...
                num2str(region_end) ' ' sample_names{j} ' 0 + ' num2str(region_start) ' ' num2str(region_start) ' 0,0,255'];
        elseif ~isempty(find(strcmp([segments_samples_arm_del_in_region(:,j); segments_samples_arm_amp_in_region(:,j)],'AMP')))
            str1=[str1 '\n' 'chr' num2str(chr) ' ' num2str(region_start) ' ' ...
                num2str(region_end) ' ' sample_names{j} ' 0 + ' num2str(region_start) ' ' num2str(region_start) ' 255,0,0'];
        end
        
        if ~isempty(segments_del_in_region)
            sample_names_j_in_cell=cell(length(segments_samples_del_in_region),1);
            sample_names_j_in_cell(:)=sample_names(j);
            idx_del=find(cellfun(@is_in_list, segments_samples_del_in_region, sample_names_j_in_cell));
            for k=1:length(idx_del)
                str1=[str1 '\n' 'chr' num2str(chr) ' ' num2str(segments_del_in_region(idx_del(k),1)) ' ' ...
                    num2str(segments_del_in_region(idx_del(k),2)) ' ' sample_names{j} ' 0 + ' num2str(segments_del_in_region(idx_del(k),1)) ...
                    ' ' num2str(segments_del_in_region(idx_del(k),1)) ' 0,0,255'];
            end
        else
            idx_del=[];
        end
        
        if ~isempty(segments_amp_in_region)
            sample_names_j_in_cell=cell(length(segments_samples_amp_in_region),1);
            sample_names_j_in_cell(:)=sample_names(j);
            idx_amp=find(cellfun(@is_in_list, segments_samples_amp_in_region, sample_names_j_in_cell));
            for k=1:length(idx_amp)
                str1=[str1 '\n' 'chr' num2str(chr) ' ' num2str(segments_amp_in_region(idx_amp(k),1)) ' ' num2str(segments_amp_in_region(idx_amp(k),2)) ' ' sample_names{j} ...
                    ' 0 + ' num2str(segments_amp_in_region(idx_amp(k),1)) ' ' num2str(segments_amp_in_region(idx_amp(k),1)) ' 255,0,0'];
            end
        else
            idx_amp=[];
        end
        
        if isempty(idx_del) && isempty(idx_amp)
            %str1=[str1 '\n' 'chr' num2str(chr) ' 0 0 ' sample_names{j} ' 0 + 0 0 0,0,0'];
            str1=[str1 '\n' 'chr' num2str(chr) ' 0 1 ' sample_names{j} ' 0 + 0 0 0,0,0'];
        end
    end
    str1=[str1 '\n' 'track name="SNPs on chip" itemRgb="On" visibility="dense"'];
    idx=find(chip_annot.chr_vec==chr & chip_annot.chr_loc_vec >= region_start & chip_annot.chr_loc_vec <= region_end);
    for j=1:length(idx)
        str1=[str1 '\n' 'chr' num2str(chr) ' ' num2str(chip_annot.chr_loc_vec(idx(j))) ' ' num2str(chip_annot.chr_loc_vec(idx(j))+1)];
    end
    file_name=fullfile(ucsc_mounted_folder,[genes_unique{i} '_' chip_annot.chip '_ucsc' normal_disease_suffix '.txt']);
    fout=fopen(file_name, 'wt');
    if fout==-1
        err=[file_name ' may already be opened.'];
        return;
    end
    fprintf(fout, '%s', sprintf(str1) );
    fprintf(fout, '\n');
    fclose(fout);
end

function flag=is_in_list(x,y)
flag=~isempty(find(strcmp(x,y)));

