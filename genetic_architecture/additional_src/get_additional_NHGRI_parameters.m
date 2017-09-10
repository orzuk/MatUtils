% Intersect with NHGRI catalog to get misssing data
% 
% Input: 
% data 
% gwas_database_name
% gwas_NHGRI_database_name
%
% Output: 
% data - gwas data structure with additional parameters 
% 
function data = get_additional_NHGRI_parameters(data, ...
    gwas_database_name, gwas_NHGRI_database_name)

num_snps = length(data.SNPs);
NHGRI_data = load(gwas_NHGRI_database_name);


[common_IDs I J] = intersect_all(empty_cell_to_empty_str(data.PUBMEDID), ...
    empty_cell_to_empty_str(NHGRI_data.PUBMEDID)); % Take common indices. Problem: it also copies empty
non_empty_inds = find(~isempty_cell(common_IDs));
I = I(non_empty_inds); J = J(non_empty_inds);
data.Study = empty_cell_to_empty_str(cell(num_snps,1));
data.Journal = data.Study;
data.Date = data.Study;
data.First_Author = data.Study;
if(~isfield(data, 'Gene'))
    data.Gene = data.Study; % For gene we already have stuff ...
end
data.Study(I) = NHGRI_data.Study(J);
data.Journal(I) = NHGRI_data.Journal(J);
data.Date(I) = NHGRI_data.Date(J);
data.First_Author(I) = NHGRI_data.First_Author(J);


field_names = fieldnames(data);
for i=1:length(field_names)
    eval(['cell_flag = iscell(data.' field_names{i} ');']); % check if cell
    if(cell_flag)
        eval_str = ['data.' field_names{i} ' = empty_cell_to_empty_str(data.' field_names{i} ');'];
        eval(eval_str);
    end
end

for i=1:num_snps
    data.SNPs{i} = str2word('-', data.SNPs{i}, 1);
    if(~isempty(data.SNPs{i}))
        if((data.SNPs{i}(end) < '0') || (data.SNPs{i}(end) > '9'))
            data.SNPs{i} = data.SNPs{i}(1:end-1);
        end
    end
end

if(~isfield(data, 'pos'))
    [data.chr data.pos] = rs_ids_to_genomic_coordinates(data.SNPs , 'hg18'); % heavy ...
    save(file_name_to_mat(gwas_database_name), '-struct', 'data');
end
for i=1:length(data.Pos) % overwrite positions when we have information
    if(~isempty(data.Pos{i}))
        tmp_pos = str2nums(data.Pos{i},1);
        if(~isempty(tmp_pos))
            data.pos(i) = tmp_pos(1);
        end
    end
    if(~isempty(data.Chr{i}))
        data.chr(i) = chr_str2num(data.Chr{i}, 'hg18');
    end
end


data.Region = data.Study;
[common_SNPs I J] = intersect_all(empty_cell_to_empty_str(data.SNPs), ...
    empty_cell_to_empty_str(NHGRI_data.SNPs)); % Take common indices
data.Region(I) = empty_cell_to_empty_str(NHGRI_data.Region(J)); % region is more complicated as it is SNP-specific!

II = setdiff(1:num_snps, I);
for i=1:length(II) % Set different regions for all other SNPs
    if(~isempty(data.Pos{II(i)}))
        data.Region{II(i)} = data.Pos{II(i)};
    else
        data.Region{II(i)} = num2str(data.pos(II(i)));
    end
end
for i=1:length(I)
    if(isempty(data.Gene{I(i)}))
        data.Gene{I(i)} = NHGRI_data.Reported_Gene_s_{J(i)}; % copy gene names from NHGRI database
    end
end

