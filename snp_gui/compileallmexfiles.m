% Compile all the mex file needed for the hmm_chrom program

cd(fullfile(matlab_libs_root_dir, 'stats/hmp')); 
mex TrainHMMFromDataEM.cpp hmm_chrom_funcs.cpp
mex FindBestPathViterbi.cpp hmm_chrom_funcs.cpp
mex SimulateSequenceFromModel.cpp hmm_chrom_funcs.cpp
mex ComputeHMMKLDistance.cpp hmm_chrom_funcs.cpp
mex ComputeHMMLogLikelihood.cpp hmm_chrom_funcs.cpp
mex ReadGenotypeFile.cpp hmm_chrom_funcs.cpp hapmap_funcs.cpp
cd(fullfile(matlab_libs_root_dir, 'stats/mog')); 
mex MixtureOfGaussiansGivenInit.cpp ../hmp/hmm_chrom_funcs.cpp
mex MixtureOfGaussiansGivenInitSingle.cpp ../hmp/hmm_chrom_funcs.cpp
