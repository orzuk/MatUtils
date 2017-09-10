function sample = sample_bnet_nocell(bnet, nsamples, varargin)
% sample_bnet_nocell Generate a random sample from a Bayes net, of size
% nsamples
% SAMPLE = SAMPLE_BNET(BNET, ...)
%
% sample(i,j) contains the value of the i'th node in the j'th sample
% Nodes are sampled in the order given by bnet.order.
% Note : This is a special-purpose sampling function for fast sampling of
% small discrete BNT. It is not good for the general case.
%

% set defauly params
n = length(bnet.dag);
sample = cell(n,1);

% get optional params
args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i},
   case 'evidence',    sample = args{i+1}(:);
   otherwise, error(['unrecognized argument ' args{i}])
  end
end

for j=bnet.order(:)'
  if isempty(sample{j})
    %ps = parents(bnet.dag, j);
    ps = bnet.parents{j};
    e = bnet.equiv_class(j);
    sample{j} = sample_node(bnet.CPD{e}, sample(ps));
  end
end
