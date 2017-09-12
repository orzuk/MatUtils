% Written by Or Zuk 5/2007
%
% Generate a toy model, simulate examples from it and then try to learn it
% from the data again and see the errors we get.
% This is done by calling the function 'DetermineModelStability'
num_points_vec = 100:200:500;
model_flag = 1; % 0 - standard HMP, 1 - spoecial SNPs model
learn_flag = 0;  % 1 : learn a new model, 0 - just infer the hidden values from true model.
learn_equ_flag = 0; % Should we learn an equivalent ('standard') model for the total copy number?
num_iters = 2; num_EM_iters = 10; num_EM_starting_points = 1;  % Parameters for EM algorithm used for learning.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set HMP Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HMM_MODEL = {};
HMM_MODEL.SPECIAL_MODELS_FLAG =1; % use SNPs special model
HMM_MODEL.x_dim = 3; HMM_MODEL.x_dim2 = 2;
%% HMM_MODEL.M = [0.999  0.001 0; 0.001 0.998 0.001; 0 0.001 0.9990]; %[0.99 0.008 0.002; 0.001 0.992 0.007; 0.009 0.001 0.99]; % Markov Model of the Mixture
%% HMM_MODEL.M = [0.99  0.005 0.005; 0.005 0.99 0.005; 0.005 0.005 0.99]; %[0.99 0.008 0.002; 0.001 0.992 0.007; 0.009 0.001 0.99]; % Markov Model of the Mixture
%% HMM_MODEL.M = [0.2  0.4 0.4; 0.4 0.2 0.4; 0.2 0.4 0.4]; %[0.99 0.008 0.002; 0.001 0.992 0.007; 0.009 0.001 0.99]; % Markov Model of the Mixture
%% HMM_MODEL.M = [0.99 0.004 0.005 0.001; 0.004 0.99 0.005 0.001; 0.001 0.004 0.99 0.005; 0.001 0.002 0.007 0.99]; % Markov Model of the Mixture
HMM_MODEL.M = [0.4  0.3 0.3; 0.3 0.4 0.3; 0.3 0.3 0.4]; %[0.99 0.008 0.002; 0.001 0.992 0.007; 0.009 0.001 0.99]; % Markov Model of the Mixture
HMM_MODEL.N =  ones(HMM_MODEL.x_dim + HMM_MODEL.SPECIAL_MODELS_FLAG*(HMM_MODEL.x_dim-1), 1); % Currently we have no mixtures - just plain Gaussians
HMM_MODEL.MU =  [0.2 0.8 1.9 2.45 3]'; %%% [-0.0462635401688274, 0.584605909385857 0.3412454325]';
HMM_MODEL.MU(3) = 2*HMM_MODEL.MU(2)-HMM_MODEL.MU(1);
HMM_MODEL.MU(4) = HMM_MODEL.MU(2)+HMM_MODEL.MU(3)-HMM_MODEL.MU(1);
HMM_MODEL.MU(5) = 2*HMM_MODEL.MU(3)-HMM_MODEL.MU(1); % Make Affine MU. Assumption: mew(k) = a+b*k
HMM_MODEL.MU = HMM_MODEL.MU(1:length(HMM_MODEL.N)); % ke it the correct dimension
HMM_MODEL.SIGMA =  0.000125*ones(1,length(HMM_MODEL.N)); % Take a Very Small St.d.
HMM_MODEL.PI = [0.2, 0.2 0.6]'; % PI has the size of x_dim (stationary on each chrom)
HMM_MODEL.LogScore = 1;  % the output score
HMM_MODEL.M_UPPERBOUNDS = [1, 0.35; 0.35, 1]; % Upper bounds on the learned transition probabilities
HMM_MODEL.USE_BOUNDS = 0; % Do not use bounds in learning
HMM_MODEL.PLACE_FLAG =0; % Use a different transition matrix for each place !!! (In special model this is true only for the LD matrices)




% Prepare the transition matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Set HMP PLace-Dependent Parameters (specifically PLACE_M matrix)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if(model_flag)
    SMALL_M = [ 0.9 0.1 ; 0.1 0.9];
    SMALL_M = [ 0.5 0.5 ; 0.5 0.5];
    %SMALL_M = [ 0.999 0.001 ; 0.001 0.999];
    HMM_MODEL.PLACE_M = zeros(num_points_vec(end), HMM_MODEL.x_dim2, HMM_MODEL.x_dim2);
    RAND_BIN_VEC = rand(1,num_points_vec(end)) < 0.5;
    for i=1:num_points_vec(end)
        HMM_MODEL.PLACE_M(i,:,:) = SMALL_M .* RAND_BIN_VEC(i) + SMALL_M([2,1],:) .* (1-RAND_BIN_VEC(i));
%%        HMM_MODEL.PLACE_M(i,1,1) = rand(1); HMM_MODEL.PLACE_M(i,1,2) = 1-HMM_MODEL.PLACE_M(i,1,1);
%%        HMM_MODEL.PLACE_M(i,2,1) = rand(1); HMM_MODEL.PLACE_M(i,2,2) = 1-HMM_MODEL.PLACE_M(i,2,1);
    end
    HMM_MODEL.PLACE_M = reshape(HMM_MODEL.PLACE_M,num_points_vec(end),HMM_MODEL.x_dim2^2); % Change to a 2-dim array
    %% TEMP = HMM_MODEL.PLACE_M(:,2); HMM_MODEL.PLACE_M(:,2) = HMM_MODEL.PLACE_M(:,3); HMM_MODEL.PLACE_M(:,3) = TEMP; %% Need to transpose? 
    HMM_MODEL.PLACE_M = HMM_MODEL.PLACE_M';  %% Is this transpose needed ? We need to check!

    % Take hapmap data !!!
    % load(  ['LD_' pop_str_vec{hapmap_population} '_' chip_type '_' genome_assembly]); % 'LD_all_chroms.mat'
    % derich = 1/120.0; % Relative Derichlet correction
    % HMM_MODEL.PLACE_M = zeros(4,length(LD{1}.Hind.PairMat{1}));
    % HMM_MODEL.PLACE_M(1,:) = (LD{1}.Hind.PairMat{1}+derich) ./ (LD{1}.Hind.PairMat{1}+LD{1}.Hind.PairMat{2}+2*derich);
    % HMM_MODEL.PLACE_M(3,:) = 1-HMM_MODEL.PLACE_M(1,:);  % Switch 2 and 3 for the C function
    % HMM_MODEL.PLACE_M(2,:) = (LD{1}.Hind.PairMat{3}+derich) ./ (LD{1}.Hind.PairMat{3}+LD{1}.Hind.PairMat{4}+2*derich);
    % HMM_MODEL.PLACE_M(4,:) = 1-HMM_MODEL.PLACE_M(2,:);
    %

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Call DetermineModelStability function which does all the job ..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ttt = cputime;
COMP_HMM_MODELS = DetermineModelStability( HMM_MODEL, num_points_vec, num_iters, ...
    num_EM_iters, num_EM_starting_points, ...
    learn_flag, learn_equ_flag);
DetStabTime = cputime - ttt


%% Simulate sequence in this function. Here we got to give it also the many transition matrices
%% [X_VEC OUT_VEC OUT_VEC_B] = SimulateSequenceFromModel(HMM_MODEL.PI, HMM_MODEL.M, HMM_MODEL.N, ...
%%     HMM_MODEL.MU, HMM_MODEL.SIGMA, HMM_MODEL.PLACE_FLAG, ...
%%     HMM_MODEL.PLACE_M, HMM_MODEL.SPECIAL_MODELS_FLAG, num_points_vec(end));
%%
%% % Transfer X_VEC to its 'real' values
%% X_VEC_alpha_genotype = bitget(X_VEC,1);
%% X_VEC_beta_genotype = bitget(X_VEC,5);
%% X_VEC_alpha_copynumber = bitget(X_VEC,2) + 2*bitget(X_VEC,3) + 4*bitget(X_VEC,4);
%% X_VEC_beta_copynumber = bitget(X_VEC,6) + 2*bitget(X_VEC,7) + 4*bitget(X_VEC,8);


%% Now plot the simulation results
%% % % figure; hold on;
%% % % plot(X_VEC_alpha_copynumber+0.01, '.');  plot(X_VEC_beta_copynumber-0.01, 'r.');
%% % % plot(OUT_VEC, 'g.');  plot(OUT_VEC_B, 'c.');
%% % % legend('alpha copy number', 'beta copy number', 'A Intensity', 'B Intensity');
%% % % xlabel('SNP location'); ylabel('Values');