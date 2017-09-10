% Intersect sets but keep all indices of elements from b mapping to a
% 
% Input: 
% a - first set
% b - second set
% 
% Output: 
% inter - intersection set
% I - indices of a in intersection 
% J - indices of b in intersection
% 
function [inter, I, J] = intersect_all(a, b) 

n = length(a); m = length(b); 
inter = []; I = []; J = [];

cover_inds = []; uncover_inds = setdiff(1:n, cover_inds);
I1 = 1;
while(~isempty(I1))
    [tmp_int, I1, J1] = intersect(a(uncover_inds), b);
    I = [I vec2row(uncover_inds(I1))]; J = [J vec2row(J1)];
    inter = [inter vec2row(a(uncover_inds(I1)))];
    uncover_inds = setdiff(1:n, I);
    
    %    b = b(setdiff(1:length(b)), J1); % remove parts already here
end
if(iscolmat(a) && iscolmat(b))
    inter = vec2column(inter);
    I = vec2column(I);
    J = vec2column(J);
end


