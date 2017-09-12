% Give a score for a pwm for how palindromic is it
Assign24MammalsGlobalConstants();

pwms_file = '../data/pwms_union.mat'; 

metric = PEARSON; % PEARSON; % EUCLIDIAN; % DOTPROD;

pal_score = pwms_palindromic_score(pwms_file, metric);

figure; hold on; hist(pal_score, 100); title('pwms palindromic scores'); 
xlabel('palindromic score'); ylabel('freq.'); 

special_pwms = {'P53', 'NRS'};
load(pwms_file);
for i=1:length(special_pwms)
	pwms_inds = strfind_cell(pwms(:,1), special_pwms{i})
	if(~isempty(pwms_inds))
		plot(pal_score(pwms_inds), zeros(length(pwms_inds),1), [color_vec(i+1) '*']);
	end
end
legend(['all-pwms' special_pwms]);
my_saveas(gcf, fullfile(html_outdir, 'pwms', 'pwms_palindromic_scores'), 'jpg'); 

for i=1:size(pwms,1)
    p_ind = strfind(pwms{i,4}, 'palindromic-score');
    if(~isempty(p_ind))
        pwms{i,4} = pwms{i,4}(1:p_ind(1)-1);
    end
    if(isempty(strfind('palindromic-score:', pwms{i,4})))
        pwms{i,4} = [pwms{i,4} '  palindromic-score: ' num2str(pal_score(i))];
    end
end
% for i=1:size(pwms,1)
%     pwms{i,4} = strdiff(pwms{i,4}, ['  palindromic-score: ' num2str(pal_score(i))]);
% end
save(pwms_file, 'pwms'); 
save_pwms_to_txt_file(pwms, [remove_suffix_from_file_name(pwms_file) '.txt'], 0);

ic_score = pwms_information_content(pwms(:,2))

pal_cutoff = 0.8;
ic_cutoff = 8 % min(ic_score)
[half_pwms half_inds half_pwms_names] = ...
    pwms_cut_in_half(pwms(1:end,:), metric, pal_cutoff, ic_cutoff)

half_pwms = [half_pwms_names half_pwms pwms(half_inds,3:4)]; pwms = half_pwms; 
save([remove_suffix_from_file_name(pwms_file) '.halfs.mat'], 'pwms'); 
save_pwms_to_txt_file(half_pwms, [remove_suffix_from_file_name(pwms_file) '.halfs.txt'], 0);


