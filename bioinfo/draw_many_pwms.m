% Display many pwms in the same plot
%
% Input: 
% pwms - vector of pwms
% counts_vec - how many sequnences generated each pwm
% w - width in figure
% l - length in figure
% pwms_names - names of pwms
% 
function Dummy = draw_many_pwms(pwms, counts_vec, w, l, pwms_names, varargin)

if(iscell(pwms))
    num_pwms = length(pwms);
else
    num_pwms = size(pwms, 3);
end
if(~exist('pwms_names', 'var'))
    for i=1:num_pwms
        pwms_names{i} = ['# seqs ' num2str(counts_vec(i))];
    end
end
for i=1:num_pwms
    if(mod(i, w*l) == 1)
        figure;
    end
    subplot(w, l, mod(i-1, w*l)+1);
    if(iscell(pwms))
        cur_pwm = pwms{i};
    else
        cur_pwm = pwms(:,:,i);
    end
    draw_pwm_new(cur_pwm, mod(i, w*l) == 1); title(pwms_names{i});
%     if(mod(i, w*l) == 1)
%         legend('A', 'C', 'G', 'T');
%     end
end

Dummy=0;



