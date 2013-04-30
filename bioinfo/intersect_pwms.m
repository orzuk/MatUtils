% Intersect two sets of pwms 
% 
% Input: 
% pwms1 - first set of pwms
% pwms2 - second set of pwms
% tolerance - how much difference do we allow between two pwms
%
% The output: 
% pwms_inter - pwms appearing in both sets
%
function pwms_inter = intersect_pwms(pwms1, pwms2, tolerance)
AssignStatsConstants;

if(~exist('tolerance', 'var')) % set default error tolerance 
    tolerance = 0.000001;
end
if(ischar(pwms1))
    if(strcmp(suffix_from_file_name(pwms1), 'mat'))
        load(pwms1);     pwms1 = pwms;
    else
       pwms1 = load_pwms_from_txt_file(pwms1);
    end
end
if(ischar(pwms2))
    if(strcmp(suffix_from_file_name(pwms2), 'mat'))
        load(pwms2);     pwms2 = pwms;
    else
       pwms2 = load_pwms_from_txt_file(pwms2);
    end
end

pwms_names1 = strrep(strrep(upper(pwms1(:,1)), '_', '.'), '.0', '.'); % First intersect names
pwms_names2 = strrep(strrep(upper(pwms2(:,1)), '_', '.'), '.0', '.');
[intersect_names inter_inds1 inter_inds2] = intersect(pwms_names1, pwms_names2);
pwms_lens1 = length_cell(pwms1(:,2)); pwms_lens2 = length_cell(pwms2(:,2));

num_shared_pwms = length(intersect_names); % Now loop over them to see that the matrices match 
pwms_dist = zeros(num_shared_pwms,1);
shared_vec = zeros(num_shared_pwms,1);
for i=1:num_shared_pwms
    i_is = i
    [pwms_dist(i) shift] = TwoPWMsToSimilarityMatrix(derich_correct(pwms1(inter_inds1(i),2), 0), ...
        derich_correct(pwms2(inter_inds2(i),2), 0), EUCLIDIAN, 1, [], 5); % metric=EUCLIDIAN, rev_comp_flag=1, min_overlap=5
    
    %    if(pwms_lens1(inter_inds1(i)) == pwms_lens2(inter_inds2(i)))
    [alpha ss] =  match_derich_correct( ...
        pwms1{inter_inds1(i),2}(:, max(1, shift+1):min(pwms_lens1(inter_inds1(i)), shift+pwms_lens2(inter_inds2(i))  )) , ...
        pwms2{inter_inds2(i),2}(:, max(1, 1-shift):min(pwms_lens2(inter_inds2(i)), -shift+pwms_lens1(inter_inds1(i))) )); % Find best alpha for a given shift
    if(i == 66)
        xxx = 999;
    end
    [pwms_dist(i) shift] = TwoPWMsToSimilarityMatrix(derich_correct(pwms1(inter_inds1(i),2), 1*alpha), ...
        derich_correct(pwms2(inter_inds2(i),2), 0*alpha), EUCLIDIAN, 1, [], 5); % metric=EUCLIDIAN, rev_comp_flag=1, min_overlap=5

    
    shared_vec(i) = pwms_dist(i) > -tolerance; % Euclidian similarity is -distance
end
shared_vec = find(shared_vec); 
pwms_inter = pwms1(inter_inds1(shared_vec),:); % perform final intersection 
figure; hist(pwms_dist, 50); 


