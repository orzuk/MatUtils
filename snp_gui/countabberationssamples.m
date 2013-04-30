% Written by Libi Hertzberg and Or Zuk 9/2007
%
% The function goes over the display file and count how many abberations were, 
% what is the total proportion of abberations, statistics on their lengths
% etc. The output of this function should be used to optimize performance
% by adjusting various parameters of the applications (say, for example we
% want 5% of the genome to be abberant - we will use this function to see
% how much is indeed abberatn in order to optimize the code of other
% functions)
%
% Output the structure DispStruct, with the following fields:
% NumCopyChanges 
% CopyNumCountsMat
function [NumCopyChanges CopyNumCountsMat] = CountAbberationsSamples(user_dir, SampleNames, chip_type, chroms)

display_dir = [user_dir  '\display'];

num_samples = length(SampleNames); 

NumCopyChanges = zeros(1,num_samples);
CopyNumCountsMat = zeros(num_samples, 5);

for cur_sample = 1:num_samples
    doing_sample = cur_sample
    sample_name = SampleNames{cur_sample}; 
    display_file_name = fullfile('display', [sample_name '_' chip_type '_disp']);
    load(fullfile(user_dir, display_file_name),'DispStruct'); % load the display structure

    [NumCopyChanges(cur_sample) CopyNumCountsMat(cur_sample,:)] = CountAbberations(DispStruct, chroms);
end

