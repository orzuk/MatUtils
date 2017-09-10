% Cut pwms to halfs and return the halfs passing some threshold
% 
% Input: 
% pwms - a set of pwms
% metric - what similarity score to use
% pal_cutoff - cutoff on palindromic score
% ic_cutoff - cutoff on information content
% 
% The output: 
% half_pwms - pwms containing halfs of some of original pwms
% 
function [half_pwms half_inds half_pwms_names] = ...
    pwms_cut_in_half(pwms, metric, pal_cutoff, ic_cutoff)
AssignGeneralConstants;

if(ischar(pwms)) % read from file
    load(pwms);
end
if(size(pwms,2) == 4) % remove annotations
    half_pwms_names = pwms(:,1); 
    pwms = pwms(:,2);
end
num_pwms = length(pwms);
pal_scores = pwms_palindromic_score(pwms, metric); 


ctr=0; half_pwms = {}; half_sides = []; half_inds = [];
for i=1:num_pwms
    ic_vec = 2-entropy(pwms{i}); % get vector
%    diff_ic_vec = diff(ic_vec); % take difference
    low_inds = find(ic_vec < 0.5); 
    high_inds = find(ic_vec > 1); 
    for j=vec2row(low_inds)
        if( ~isempty(high_inds) && (j < max(high_inds)) && (j > min(high_inds)) ) % found a 'dip' in ic vec
            cut_ind = j;
            half_inds = [half_inds i i];
            half_pwms{ctr+1} = pwms{i}(:,1:cut_ind); half_sides(ctr+1) = LEFT;
            half_pwms{ctr+2} = pwms{i}(:,cut_ind:end); half_sides(ctr+1) = RIGHT;
            ctr=ctr+2;
            break;
        end
    end
end
    
% ic_scores = pwms_information_content(pwms); 

pal_inds = find(pal_scores > pal_cutoff);
pal_inds2 = vec2row(reshape(repmat(pal_inds, 1, 2)', size(pal_inds,1)*2,1));
half_inds = [vec2row(half_inds) vec2row(pal_inds2)];
half_sides = [vec2row(half_sides) vec2row(pal_inds2)];
half_pwms = [half_pwms cell(1, length(pal_inds2))];
for i=1:length(pal_inds)
    L = size(pwms{half_inds(i)}, 2); % get pwms length
    half_pwms{ctr+2*i-1} = pwms{half_inds(i)}(:,1:ceil(L/2)); % take left part  
    half_sides(ctr+2*i-1) = LEFT;
    half_pwms{ctr+2*i} = pwms{half_inds(i)}(:,floor(L/2):end); % take right part  
    half_sides(ctr+2*i) = RIGHT;
end

half_ic_scores = pwms_information_content(half_pwms);
half_ic_inds = find(half_ic_scores > ic_cutoff);
half_pwms = half_pwms(half_ic_inds); half_inds = half_inds(half_ic_inds);
half_sides = num2direction(half_sides(half_ic_inds));
half_pwms_names = half_pwms_names(half_inds);
for i=1:length(half_pwms_names)
    half_pwms_names{i} = [half_pwms_names{i} '.' half_sides{i}];
end
    
% Get rid of matrices appearing twice (once from palindromic and once from 
% information dip)
[half_pwms_names u] = unique(half_pwms_names);
half_inds = half_inds(u); 
half_pwms = vec2column(half_pwms(u)); 


