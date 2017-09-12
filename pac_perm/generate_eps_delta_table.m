% Calculate the values of epsilon and delta and number samples we want for 
% various different numbers. 

% Choose if to take both top and bottom genes or just top
global ONE_SIDE TWO_SIDES;
ONE_SIDE = 1; TWO_SIDES = 0;
one_side_flag = TWO_SIDES;

% Choose if to compute overlap between true and sampled or between two
% sampled
global TRUE_AND_SAMPLED TWO_SAMPLED;
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
true_corr_flag  = TWO_SAMPLED;


alpha = 0.0046; % 70 genes in VantVeer new 
%%%alpha = 0.012;  % 70 genes in VantVeer old

%%%%alpha_vec =  [ 0.012 0.00403001 0.012]; % corresponds to ~70 and ~700 genes


global UNIFORM GAUSSIAN LINEAR FROM_DATA mix_GAUSSIAN student_t;
UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2; FROM_DATA = 3; mix_GAUSSIAN=4;
student_t=5;% in the last one we take q simply according to the data bins
rand_flag = GAUSSIAN;

% data files
OLD_VANT_VEER = 1; NEW_ROSSETA = 2; WANG = 3;     % Breast gene expression
TOPIC = 4; % Author-topic matching
NIPS_DOROTHEA = 5; NIPS_ARCENE = 6; NIPS_DEXTER = 7; NIPS_GISETTE = 8; % Various datasetes from NIPS contest
ROSEN = 9; HEPATO =10; % New data's from Assif
Bhattacharjee=11;  BEER=12; % New Lung data
YEOH = 13; % leukemia
BRAIN = 14; KIM = 15; % New aging related datasets

RAND_DATA = 16; % Here we randomize a data so that we have control on it, and also can work with matlab 6.5 (no loading from files)


load('all_datas_vars.mat'); % Load the data struct

data_flag =  WANG; % OLD_VANT_VEER; % NEW_ROSSETA; % NEW_ROSSETA; % NEW_ROSSETA;  % Choose which data to plot


% Here we have at least 1-eps overlap, with probability 1-delta
 eps_vec = [0.98 0.95 0.9 0.8 0.5]; delta_vec = [0.1 0.5]; 
%eps_vec = [0.9998]; delta_vec = [0.2]; 

nsamples_matrix = zeros(length(eps_vec), length(delta_vec)); % Matrix storing the number of samples needed

Ngenes = data_output_strct{data_flag}.Ngenes; % This should be taken from the data, after preprocessing 
sigma = sqrt(data_output_strct{data_flag}.Var_C); % Here determine the st.d. at infinity. Try taking the sqrt  
miu=0; prior = []; 

for i=1:length(eps_vec)
    for j=1:length(delta_vec)
        nsamples_matrix(i,j) = ...
            compute_num_samples_needed(rand_flag, one_side_flag, true_corr_flag, sigma, alpha, eps_vec(i), delta_vec(j), Ngenes, miu, prior)
    end
end
        

