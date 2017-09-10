% This function finds which of the samples are considered 'Normals'.
% The convention is that the last two characters in a Normal sample 
% must be '_n' for it to be considered 'Normal'. Otherwise it is not
% considered 'Normal'
%
% Input: samples - the names of the sample
%
% Output: normal_ind - Indexes of the normal samples
function normal_ind = get_normal_copy_num_samples_ind(samples)

num_samples = length(samples);
normal_ind = zeros(1, num_samples);
for i = 1:num_samples
    sample = char(samples{i});
    underscore_ind = strfind(sample, '_');
    if(length(underscore_ind))
        if(strcmp('_n', sample(underscore_ind(end):underscore_ind(end)+1)))
            normal_ind(i) = 1;
        end
    end
end

normal_ind = find(normal_ind==1);