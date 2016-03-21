% Read Variant Call Format file to .mat
%
% Input:
% vcf_file_name - input file name
% output_file - where to save .mat version
%
% Output:
% VCF - structure with all relevant data
%
function VCF = read_vcf(vcf_file_name, output_file)

skip_lines = '##'; % 10
[VCF, R, R_skipped] = ReadDataFile(vcf_file_name, output_file, -1, skip_lines); % skip first 10 lines

ctr=1;
for i=1:size(R_skipped, 1)
    if(~isempty(strfind(R_skipped{i,1}, 'INFO')))
        ID_str = strsplit(R_skipped{i,1}, '=');
        VCF.field_names{ctr} = str2word(',', ID_str{3}, 1);
        VCF.field_descriptions{ctr} = strdiff(strdiff(strdiff(ID_str{end}, '"'), ''''), '>');
        ctr=ctr+1;
    end
end

VCF.sample_names = R(1,10:end);
R = R(2:end,10:end); % get only relevant indices 
num_people = size(R,2);
num_snps = size(R,1);

VCF.GenotypeMat = zeros(num_snps, num_people)-1; % minus one indicates no genotype
VCF.CoverageMat = zeros(num_snps, num_people);
VCF.Population = cell(num_people,1);

for i=1:num_people
    VCF.Population{i} = str2word('.', remove_dir_from_file_name(vcf_file_name), 1);
end
for i=1:num_snps
    if(mod(i,100) == 0)
        read_snp_in_VCF = i
    end
    for j=1:num_people
        %        j_is = j
        tmp_nums = str2nums(R{i,j});
        if(length(tmp_nums) == 3)
            VCF.GenotypeMat(i,j) = tmp_nums(1) + tmp_nums(2);
        end
        if(~isempty(tmp_nums))
            VCF.CoverageMat(i,j) = tmp_nums(end);
        end
    end
end
if(isempty(output_file))
    output_file = file_name_to_mat(vcf_file_name);
end
save(output_file, '-struct', 'VCF');
