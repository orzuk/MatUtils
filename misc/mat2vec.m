% Convert matrix to a vector (written by Libi Hertzberg)
function vec = mat2vec(mat)

%vec = mat(:);

vec = mat'; vec = vec(:); % NEW! make it compatible with matlab's vec2mat function ! 