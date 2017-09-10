% Parse simulation data from sfs-code program of Hernansez et al. 
%
% Input: 
% sfs_output_file_name - name of sfs output file 
%
% Output: 
% sfs_struct - structure with data representing alleles, frequencies etc. 
% 
function [sfs_struct] = parse_sfs_code_output ( sfs_output_file_name )


if(~exist('sfs_output_file_name', 'var') || isempty(sfs_output_file_name))
    sfs_output_file_name = 'C:/users/user/Downloads/sfscode/bin/out';
end

% write parser here 
R = loadcellfile(sfs_output_file_name, [], ';'); % read cell file

command_str = R{1,1}; 
seed = str2nums(R{2,1});
iterations = str2nums(R{3,1}, [], [], 0); % get number of iterations (could be a vector !!!) 
N_final = str2nums(R{4,1}, [], [], 0);
Males = str2nums(R{5,1}, [], [], 0)
freq_vec = zeros(size(R,2), 1); 
for i=1:size(R,2) % read frequencies 
    tmp_freq  = str2num(str2word(',', R{6,i}, 12))
    if(~isempty(tmp_freq))
        freq_vec(i) = tmp_freq; 
    end
end
freq_vec = sort(freq_vec(freq_vec > 0)) ./ (2*N_final(1)); % normalize by N_final (which one to take??) 

sfs_struct = var2struct(command_str, seed, iterations, N_final, Males, freq_vec)

