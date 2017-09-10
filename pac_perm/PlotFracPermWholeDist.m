% A script plotting the overlap distribution 
TOL = 0.000001;

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


% Flag saying if we take alpha to be const or take alpha*N ~ n
const_alpha_flag = 0; 

% Here do only maple integral to show a pic.
res = 0.005; max_sigma = 4; alpha_res = 0.01;
%alpha_vec = [alpha_res:alpha_res:0.5];
sigma_vec = [res:res:max_sigma];

dist_std = 1; % The std. of the original distribution !!! 

%JointDensFrac(y, sig, C_alpha, one_side_flag)


% Here there's only one alpha!!!
alpha = 0.012; N_g = 10000;

%%%%% Maple Integral 
%%%% int( (1/sqrt(2*Pi)) *  exp(-c*c/2) * int( (1/(sigma*sqrt(2*Pi))) * exp(-x*x/(2*sigma*sigma)),x=x_alpha-c..infinity), c=c_alpha..infinity);
C_alpha = norminv(1-0.5*alpha);  % Get the C_alpha vector


miu = []; prior = [];

f_mean_vec = zeros(1,length(sigma_vec));
f_std_vec = zeros(1,length(sigma_vec));

f_vec = [res:res:1-res];

ttt = cputime;

    nsamples_vec = 1./sigma_vec.^2 + 3;

length(sigma_vec)
for i = 1:length(sigma_vec)
    % Choose kind of scaling to take:
    if(const_alpha_flag)
        cur_alpha=alpha;
    else
        cur_alpha = min(nsamples_vec(i)/N_g,1.0-TOL);  % Problems if this is bigger than one!!!!!!
    end
        
    [f_mean_vec(i), f_std_vec(i), x_alpha] = compute_inf_limit_fraction(rand_flag, one_side_flag, true_corr_flag, ...
        dist_std, nsamples_vec(i), cur_alpha,miu,prior);
    do_i = i
end



% Now compute the probability distribution:
prob_f_mat = zeros(length(f_vec),length(sigma_vec));


for i=1:length(sigma_vec)
    prob_f_mat(:,i) = (sqrt(N_g)/(sqrt(2*pi) * f_std_vec(i))) .* exp (  - (f_vec - f_mean_vec(i)).^2 .* N_g ./ (2 .* f_std_vec(i).^2)  );
    do_i_2nd_loop = i
end



time_elapsed = cputime - ttt

figure; hold on; imagesc(f_vec, sigma_vec,  min(prob_f_mat',15)); colorbar;
ylabel('\Sigma_n'); xlabel('f'); AXIS([res 1-res 0 sigma_vec(end)]);
if(const_alpha_flag)
    title(['Prob. density of f for \alpha = ' num2str(alpha) ' and N_g = ' num2str(N_g)]);
else
    title(['Prob. density of f for \alpha = N_g/n and N_g = ' num2str(N_g)]);
end

