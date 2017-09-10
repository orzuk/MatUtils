% Written by Or Zuk 6/2007
%
% This function is designed to cheat the SNP_GUI and make him think that
% HMM was run for all the samples, by changing the file after_run_hmm_ ...
%
function OutputFilesNames = ...
    FakeAfterRunHMMFile(user_dir, SampleNames, LDStruct, SNPChipAnnotStruct, HMMParamsStruct)

% Various data types
HMMParamsStruct.use_affy_genotypes = 0; % Should be eliminated - input as a button ...
chip_type = lower(SNPChipAnnotStruct.chip);
num_samples = length(SampleNames);

for cur_sample = 1:num_samples % Outer loop on samples ....
    sample_name = SampleNames{cur_sample};

    display_file_name = fullfile('display',[sample_name '_' chip_type '_disp']);
    hmm_out_file_name = fullfile('hmm_out',[sample_name '_' chip_type '_hmm_out']);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set the output file names:
    OutputFilesNames.disp_files{cur_sample+1} = [display_file_name '.mat'];
    OutputFilesNames.hmm_out_files{cur_sample+1} = [hmm_out_file_name '.mat'];
    

end % loop on samples
OutputFilesNames.disp_files{1} = fullfile('display', 'AllSamplesAverage_disp.mat'); % Add also files names fore averages
OutputFilesNames.hmm_out_files{1} = fullfile('hmm_out', 'AllSamplesAverage.mat');

samples = SampleNames;


% save and run over 
save(fullfile(user_dir, ['after_run_hmm_' chip_type '.mat']), 'samples', 'OutputFilesNames');






