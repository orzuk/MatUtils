%% An interface function for calling the C function which does all the work.
%% The stationaty distribution PI is fixed after calling the C function. 
%% The function has many many input AND output variables. Below is their description:
%%
%% Output (this is the easier part ...):
%% -------------------------------------
%% PI - The stationary Markov distribution
%% M - The Markov transition matrix
%% N - The Emission/Noise matrix
%% MEW - The means of the contunuous (Gaussian) outputs
%% SIGMA - The st.ds. of the continuous (Gaussian) outputs.
%% LogScore - The log-likelihood of the optimal model
%% 
%% Input (Many many variables and flags ...)
%% -------------------------------------
%% DataVec - First data intensity (A chrom) 
%% DataVec2 - Secont data intensity (B chrom)
%% locations - Chromosomal locations
%% use_locations - Flag saying if the model should use these chromosomal locations
%% HMM_x_dim - Dimension of the Markov state space
%% HMM_y_dim - Dimension of the observations OR number of Gaussians in the observation mixture
%% place_flag - Whether to use place-dependent matrix information or not.
%% do_fold_change - should we do fold change? when? why? ???
%% mean_vec_rep - The means of the Y's if using place flag - Given as input and not learned (in our SNP case irrelevant) 
%% std_vec_rep - The st.d.s of the Y's if using place flag - Given as input and not learned (in our SNP case irrelevant)
%% PLACE_M - the place-dependent Markov matrices
%% SPECIAL_MODELS_FLAG - flag indicating if we have some taylored specific Markov Model for SNPs
%% M_UPPERBOUNDS - An upperbound on the transition matrix. This is used in order to avoid 'Singular' models 
%% USE_BOUNDS - Whether or not to use the upperbounds and avoid singularity
%% num_EM_iters - Number of EM iterations to perform
%% num_EM_starting_points - Number of different EM (random) sarting points
%% EM_tolerance - The improvement in score we require. If it's lower then we stop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [PI, M, N, MEW, SIGMA, LogScore] = TrainHMMFromDataEMMatlab(DataVec, DataVec2, locations, use_locations, ...
    HMM_x_dim, HMM_y_dim, ...
    place_flag, do_fold_change, mean_vec_rep, std_vec_rep, PLACE_M, ...
    SPECIAL_MODELS_FLAG, M_UPPERBOUNDS, USE_BOUNDS, ...
    num_EM_iters, num_EM_starting_points, EM_tolerance);       


% Fund the starts of the chromosomes
AllStarts = [1 1+find(locations(1:end-1) > locations(2:end) )];


% Calling the C mex-function 
[PI, M, N, MEW, SIGMA, LogScore] = TrainHMMFromDataEM(DataVec, DataVec2, locations, use_locations, AllStarts, ...
    HMM_x_dim, HMM_y_dim, ...
    place_flag, do_fold_change, mean_vec_rep, std_vec_rep, PLACE_M, ...
    SPECIAL_MODELS_FLAG, M_UPPERBOUNDS, USE_BOUNDS, ...
    num_EM_iters, num_EM_starting_points, EM_tolerance);    

% Now correct PI - This is no good in the C code, so we change PI to be stationary here
[eigenvecs eigenvals] = eig(M');
[val ind] = min(abs(diag(eigenvals-eye(HMM_x_dim)))); % Find the index for which the eigenvalue is closest to one : 
PI = eigenvecs(:,ind) / sum(eigenvecs(:,ind));

eps = 0.000000001; % Check that PI is summed to one
if(sum(PI) > 1 + eps) 
    ERROR_PIPIPIPIPIPI = PI
end



