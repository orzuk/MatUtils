% Flip MAF and effect size for an allele
%
% Input:
% f_vec - RAF
% grr_vec - effect size
% trait_type - binary (default) or QTL
% flip_mode - 0: flip all (default), 1: only RAF > 0.5, 2: only negative effect sizes
%
% Output:
% flipped_f_vec - flipped RAF
% flipped_grr_vec - flipped effect size
% flipped_flag_vec - binary vector indicating which were flipped (ones)
%
function [flipped_f_vec flipped_grr_vec flipped_flag_vec] = ...
    flip_allele(f_vec, grr_vec, trait_type, flip_mode)

if(~exist('trait_type', 'var') || isempty(trait_type)) % default trait type is binary
    trait_type = 'binary';
end
if(~exist('flip_mode', 'var') || isempty(flip_mode))
    flip_mode = 0;  % default: set to zero
end
switch flip_mode
    case 0 % flip ALL of them
        flipped_inds = 1:length(f_vec);        
    case 1 % flip all with RAF > 0.5
        flipped_inds = find(f_vec > 0.5);
    case 2 % flip all with effect size negative-beta or grr<1
        switch trait_type
            case {'binary', 'Binary'}
                flipped_inds = find(grr_vec < 1);
            case {'QTL', 'Quantitative'}
                flipped_inds = find(grr_vec < 0);
        end
end
flipped_f_vec = f_vec;  flipped_grr_vec = grr_vec; % start with copying input
flipped_f_vec(flipped_inds) = 1-f_vec(flipped_inds);
switch trait_type
    case {'binary', 'Binary'}
        flipped_grr_vec(flipped_inds) = 1 ./ grr_vec(flipped_inds);
    case {'QTL', 'Quantitative'}
        flipped_grr_vec(flipped_inds) = -grr_vec(flipped_inds);
end
flipped_flag_vec = zeros(size(f_vec)); flipped_flag_vec(flipped_inds) = 1; % record which ones were flipped
