function scores = my_BIC_score_dags( data_counts, nsamples, dags)
% Give the MDL-BIC  score for the dags, based on the data counts. Works
% only for small binary bnets.

num_datas = size(data_counts,1);
m=size(data_counts,2);


scores = zeros(num_datas, length(dags));

% size(data_counts)
for j=1:length(dags)
    size(data_counts(:,j))
    P = project_on_bnet(data_counts', dags{j});
%     size(P)
%     
%     PPP = P

    %     size(data_counts')
    %     size(repmat(P, 1, num_datas))
    %     size(sum(data_counts' .* log(repmat(P, 1, num_datas)) ))
    
%     sss = size(scores(:j))
%     zzz = size( sum(data_counts' .* log(P)))
    scores(:,j) = sum(data_counts' .* log(P));

%     sss = scores(:,j)
    
    scores(j,:) = scores(j,:) - dim_dag(dags{j}) * log(nsamples) / 2;


end
