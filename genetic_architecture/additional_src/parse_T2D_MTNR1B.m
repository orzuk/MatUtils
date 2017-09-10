% Analyze the rare-alleles appearing in the paper on MTNR1B rare-variants on T2D

%function parse_T2D_MTNR1B()

MTNR1B_rare_file = '../../common_disease_model/data/rare/MTNR1B_Type2Diabetes/forty_variants_in_MTNR1B.txt';

R = ReadDataFile(MTNR1B_rare_file, [], [], 5, 9); 
R.Variant = strrep_cell(R.Variant, 'NUM_PARTICIPANTS', 'Total:'); 

%function [S R] = ReadDataFile(data_file, output_file, cell_to_mat, skip_lines, delimiter, varargin)


num_cases = R.T2D(end)
num_ifg_cases = R.IFG(end)
num_controls = R.NG(end)

R.num_variants_cases = round(R.T2D(5:end-1)*num_cases*2./100) % account for two chromosomes
R.num_variants_ifg = round(R.IFG(5:end-1)*num_ifg_cases*2./100) % account for two chromosomes
R.num_variants_controls = round(R.NG(5:end-1)*num_controls*2./100) % account for two chromosomes

total_num_variants_cases = sum(R.num_variants_cases)
total_num_variants_ifg_cases = sum(R.num_variants_ifg)
total_num_variants_controls = sum(R.num_variants_controls)

total_num_variants_cases  / num_cases
total_num_variants_controls  / num_controls

(total_num_variants_cases  / num_cases) / (total_num_variants_controls  / num_controls)


LOF_inds = find(R.Loss_of_Function)-4;
total_num_LOF_cases = sum(R.num_variants_cases(LOF_inds))
total_num_LOF_controls = sum(R.num_variants_controls(LOF_inds))
total_num_LOF_ifg = sum(R.num_variants_ifg(LOF_inds))

(total_num_LOF_cases  / num_cases) / (total_num_LOF_controls  / num_controls)



R.MAF_T2D = R.T2D;
R.MAF_IFG = R.IFG;
R.MAF_NG = R.NG;


R.carriers_T2D = [ 0 0 0 0 R.num_variants_cases' 0]';
R.carriers_IFG = [ 0 0 0 0 R.num_variants_ifg' 0]';
R.carriers_NG = [ 0 0 0 0 R.num_variants_controls' 0]';
R.LossOfFunction = R.Loss_of_Function; R.MAF_controls = R.MAF_controls____;
R.Common = [1 1 1 1 zeros(1,37)]';

R = rmfield(R, {'num_variants_cases', 'num_variants_ifg', 'num_variants_controls', 'Loss_of_Function', 'MAF_controls____', ...
    'T2D', 'IFG', 'NG'});


R.carriers_T2D(end) = sum(R.carriers_T2D(1:end-1) .* (1-R.Common(1:end-1)));
R.carriers_IFG(end) = sum(R.carriers_IFG(1:end-1) .* (1-R.Common(1:end-1)));
R.carriers_NG(end) = sum(R.carriers_NG(1:end-1) .* (1-R.Common(1:end-1)));
R.carriers_T2D(end+1) = R.MAF_T2D(end);
R.carriers_IFG(end+1) = R.MAF_IFG(end);
R.carriers_NG(end+1) = R.MAF_NG(end);
R.MAF_T2D(end) = sum(R.MAF_T2D(1:end-1) .* (1-R.Common(1:end-1)));
R.MAF_IFG(end) = sum(R.MAF_IFG(1:end-1) .* (1-R.Common(1:end-1)));
R.MAF_NG(end) = sum(R.MAF_NG(1:end-1) .* (1-R.Common(1:end-1)));
R.MAF_T2D(end+1) = 0; R.MAF_IFG(end+1) = 0; R.MAF_NG(end+1) = 0;
R.Variant{end+1} = 'Individuals';
R.LossOfFunction(end+1) = 0;
R.Common(end+1)=0;
R.MAF_controls{end+1} = '';



R = orderfields(R, {'Variant', 'MAF_controls', 'MAF_T2D', 'carriers_T2D', 'MAF_IFG', 'carriers_IFG', ...
    'MAF_NG', 'carriers_NG', 'LossOfFunction', 'Common'});
WriteDataFile(R, [remove_suffix_from_file_name(MTNR1B_rare_file) '_preprocessed.txt'], 0); 




