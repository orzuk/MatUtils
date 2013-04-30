% Display first principal components of data
function [coefs,score,latent] = DisplayPCA(y, m, trans_flag, samp_n, varargin)

if(ischar(y)) % load a file
    file_name = y;
    y = load(y); 
    if(isa(y, 'struct'))
        y = y.data;
    end
% Save in matlab format 
%     if(~strcmp(suffix_from_file_name(file_name), 'mat'))
%         data = y;
%         save([remove_suffix_from_file_name(file_name) '.mat'], 'data');
%         clear data;
%     end
end

if(exist('trans_flag', 'var') && (~isempty(trans_flag)))
    if(trans_flag)
        y = y';
    end
end
if(exist('samp_n', 'var')) % sample lines
    if(samp_n < size(y,2))
        p = rand_nchoosek(size(y,2), samp_n);
        y = y(:,find(p));
    end
end
[coefs,score,latent] = princomp(y); % score is the projection of the X vectors to the PCA axes coefs

if(m == 2)
    figure; plot(score(:,1), score(:,2), '*', 'linewidth', 2);
    xlabel('1st pca'); ylabel('2nd pca'); 
else
    figure; plot3(score(:,1), score(:,2), score(:,3), '*', 'linewidth', 2);
    xlabel('1st pca'); ylabel('2nd pca'); zlabel('3rd pca'); 
end

