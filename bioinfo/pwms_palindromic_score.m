% Give a score for a pwm for how palindromic is it
% 
% Input: 
% pwms - a set of pwms
% metric - what similarity score to use
% 
% The output: 
% pal_score - score of best palindromic match
% 
function pal_score = pwms_palindromic_score(pwms, metric)

if(ischar(pwms)) % read from file
    load(pwms);
end
if(size(pwms,2) == 4) % remove annotations
    pwms = pwms(:,2);
end


num_pwms = length(pwms);
pal_score = zeros(num_pwms,1);
for i=1:num_pwms
    pwm_rev_comp = pwmrcomplement(pwms{i});
    pal_score(i) = max(TwoPWMsToSimilarityMatrix(pwms(i), {pwm_rev_comp}, metric, 0,0)); 	% pal_score(i)
end

