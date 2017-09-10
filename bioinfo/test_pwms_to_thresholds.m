% Test the function pwms_to_thresholds
load('../data/pwms_union.mat'); 
p = pwms(1:10,2); alpha = 0.001;

% for i=1:length(p)
%     p{i} = single(p{i});
% end
t = pwms_to_thresholds(p, alpha)


[t alpha_IC] = pwms_to_thresholds(p) % here let the pwms themselves determine their thresholds

genomic_model = [0.3 0.2 0.2 0.3]; % take genomic (AT-rich) background model
uniform_model = [0.25 0.25 0.25 0.25]; % take uniform background model

alpha = 0.001;
iters = 50000;
% t_genomic = pwms_to_thresholds(p, alpha, iters, genomic_model)
t_uniform = pwms_to_thresholds(p, alpha, iters, uniform_model)'
t = pwms_to_thresholds(p, alpha, iters)'
t2 = pwms_to_thresholds(p, alpha, iters)'
t_and_uniform_diff = (t - t_uniform)
t_diff = t-t2
