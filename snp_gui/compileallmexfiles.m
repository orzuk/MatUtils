% Compile all the mex file needed for the hmm_chrom program

mex TrainHMMFromDataEM.cpp hmm_chrom_funcs.cpp
mex FindBestPathViterbi.cpp hmm_chrom_funcs.cpp
mex SimulateSequenceFromModel.cpp hmm_chrom_funcs.cpp
mex ComputeHMMKLDistance.cpp hmm_chrom_funcs.cpp
mex ComputeHMMLogLikelihood.cpp hmm_chrom_funcs.cpp
mex ReadGenotypeFile.cpp hmm_chrom_funcs.cpp hapmap_funcs.cpp
mex MixtureOfGaussiansGivenInit.cpp hmm_chrom_funcs.cpp
mex MixtureOfGaussiansGivenInitSingle.cpp hmm_chrom_funcs.cpp
