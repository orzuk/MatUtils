% Choose the probability distribution of the TRUE corrleations
UNIFORM = 0; GAUSSIAN = 1; LINEAR = 2; FROM_DATA = 3;mix_GAUSSIAN=4;
student_t=5;% in the last one we take q simply according to the data bins
rand_flag = GAUSSIAN;

% data files
OLD_VANT_VEER = 1; NEW_ROSSETA = 2; WANG = 3;     % gene expression
TOPIC = 4; % Author-topic matching
NIPS_DOROTHEA = 5; NIPS_ARCENE = 6; NIPS_DEXTER = 7; NIPS_GISETTE = 8; % Various datasetes from NIPS contest
RAND_DATA = 9; % Here we randomize a data so that we have control on it, and also can work with matlab 6.5 (no loading from files)


% Vector saying which datasets come in a sparse form
IS_SPARSE_VEC = [0,0,0,0,0,0,0,0,0];
CALC_CORRS_VEC = [1,1,1,0,0,0,0,1,1];  % calc for Gizette [1,1,1,0,0,0,0,0];
CALC_VAR_VEC = [0,0,0,0,0,0,0,0,0];


% Choose if to take both top and bottom genes or just top
ONE_SIDE = 1; TWO_SIDES = 0;
one_side_flag = TWO_SIDES;

% Choose if to compute overlap between true and sampled or between two
% sampled
TRUE_AND_SAMPLED = 1; TWO_SAMPLED = 0;
true_corr_flag  = TRUE_AND_SAMPLED;


% Here do only maple integral to show a pic.
res = 0.01; max_sigma = 4; alpha_res = 0.01;
alpha_vec = [alpha_res:alpha_res:0.5];
sigma_vec = [res:res:max_sigma];

dist_std = 1; % The std. of the original distribution !!! 

%JointDensFrac(y, sig, C_alpha, one_side_flag)

%%%%% Maple Integral 
%%%% int( (1/sqrt(2*Pi)) *  exp(-c*c/2) * int( (1/(sigma*sqrt(2*Pi))) * exp(-x*x/(2*sigma*sigma)),x=x_alpha-c..infinity), c=c_alpha..infinity);
C_alpha = norminv(1-0.5*alpha_vec);  % Get the C_alpha vector


miu = []; prior = [];

f_mean_mat = zeros(length(alpha_vec),length(sigma_vec));
f_std_mat = zeros(length(alpha_vec),length(sigma_vec));

ttt = cputime;


length(sigma_vec)
for i = 1:length(sigma_vec)
    i
    for j=1:length(alpha_vec)
        %sigma_vec(i)
        %C_alpha(j)
        %            F = inline('(1/sqrt(2.*pi)) .* exp(-c.*c./2) .* (1 - normcdf( (C_alpha(j).*(1+sigma)-c)./sigma_vec(i) ))');
%%%        f_mean_mat(j,i) = (1.0/alpha_vec(j))*quadl('JointDensFrac', C_alpha(j), 99999, [], [], sigma_vec(i), C_alpha(j), one_side_flag);
        
         nsamples = 3+1/sigma_vec(i)^2;
         [f_mean_mat(j,i), f_std_mat(j,i), x_alpha] = compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, dist_std, nsamples, alpha_vec(j),miu,prior);
%         
%         
%         x_alpha = C_alpha * sqrt(1+sigma.^2);
%         f_mean_mat(j,i) = (1.0/alpha)*quadl('JointDensFrac', C_alpha, 100*(1+sqrt(sigma)), TOL, [],...
%                 sigma, C_alpha, one_side_flag);
        
    end
end

time_elapsed = cputime - ttt

figure; hold on; imagesc(alpha_vec, sigma_vec,  f_mean_mat'); colorbar;  
ylabel('\Sigma_n'); xlabel('frac. \alpha'); AXIS([alpha_res 0.5 0 sigma_vec(end)]);
title('Mean f');

figure; hold on; imagesc(alpha_vec, sigma_vec,  f_std_mat'); colorbar;  
ylabel('\Sigma_n'); xlabel('frac. \alpha'); AXIS([alpha_res 0.5 0 sigma_vec(end)]);
title('Std. f');





