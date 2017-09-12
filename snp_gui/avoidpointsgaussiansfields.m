% For each data point give a class according to which gaussian does it
% belong out of the multi-dimensional MOG.
function C = ClassifyMixtureOfGaussiansMultiDim(x, Mew_vecs, Sigma_mats, P_vec) 

num_gaussians = length(Sigma_mats);
dim=size(x,1) % Data dimension
N=size(x,2) % # data points

p = zeros(N, num_gaussians);
for m=1:num_gaussians
   p(:,m)=exp(- sum(((x'-repmat(Mew_vecs(m,:),N,1)) * inv(Sigma_mats{m})  .* ...
       (x'-repmat(Mew_vecs(m,:),N,1)   ))  ./ 2,2 ) ) ./ (sqrt((2*pi)^dim.*det(Sigma_mats{m})));
end

[DummyVal C] = max(p');