% Compute bootstrap of samples 
function [mean_corr,var_corr,f,var_f,num_of_appearances,num_of_double_appearances,hist_var] = ...
tob_removed_with_bootsrapping(data,samples,group_sizes)

cd E:\Gene_Corr\Or

alpha=0.012;

% data=gene_expression_log2;
% samples=Prognosis;
% N_TOP=267;
group_sizes=min(group_sizes):1:max(group_sizes);
N_genes=size(data,1);
N_samples=size(data,2);
N_TOP=floor(alpha*N_genes);
idx_poor=find(samples==1);
idx_good=find(samples==0);
N_poor=length(idx_poor);
N_good=N_samples-N_poor;
ratio=N_poor/N_samples;
N_iterations=200;
% min_i=ceil(1/ratio);
% max_i=floor(N_samples/2);
num_of_sizes=length(group_sizes);
mean_corr=zeros(N_genes,num_of_sizes);
var_corr=zeros(N_genes,num_of_sizes);
f=zeros(1,num_of_sizes);
var_f=zeros(1,num_of_sizes);
num_of_appearances=zeros(N_genes,num_of_sizes);
num_of_double_appearances=zeros(N_genes,num_of_sizes);
CORR=[];
for i=1:num_of_sizes
    group_size=group_sizes(i);
%     N_TOP=group_size;
    idx_groups=round(rand(N_iterations,group_size-2)*(N_samples-1))+1;
    idx_groups=[idx_groups idx_poor(round(rand(N_iterations,1)*(N_poor-1))+1)...
        idx_good(round(rand(N_iterations,1)*(N_good-1))+1)]; 
    ff=[];
    var_ff=[];
    corr=[];
    for k=1:N_iterations

        norm_data=dna_N(data(:,idx_groups(k,:)));
        curr_samples=samples(idx_groups(k,:));
        norm_samples=curr_samples-mean(curr_samples);
        norm_samples=norm_samples/sqrt(sum(norm_samples.^2));
        C=norm_data*norm_samples;
        C(C==1)=0.9999999;
        C(C==-1)=-0.9999999;
        corr(:,k)=atanh(C);
%         idx1=find(isnan(corr(:,k))~=0);
%         idx2=find(abs(corr(:,k))<1.e14);
%         if(~isempty(idx))
%             corr(idx,k);
%             norm_data(idx,:);
%             norm_samples;
%             corr(idx,k)=1000;
%         end
        CORR(end+1)=corr(100,k);
        [val,ind]=sort(-abs(corr(:,k)));
        new_list=ind(1:N_TOP);
        idx=ind(1:N_TOP);
        num_of_appearances(idx,i)=num_of_appearances(idx,i)+1;
        if(k>1)
            idx_double=intersect(old_list,new_list);
            num_of_double_appearances(idx_double,i)=num_of_double_appearances(idx_double,i)+1;
            ff(end+1)=length(idx_double);
        end
    old_list=new_list;
    end
    hist_var(i)=mean(var(corr,0,1))
    mean_corr(:,i)=mean(corr,2);
    var_corr(:,i)=mean(corr.^2,2)-mean(corr,2).^2;
    num_of_appearances(:,i)=num_of_appearances(:,i)/N_iterations;
    num_of_double_appearances(:,i)=num_of_double_appearances(:,i)/length(ff);
    f(i)=mean(ff);
    var_f(i)=mean(ff.^2)-f(i)^2;
end
f=f/N_TOP;    
var_f=sqrt(var_f)/N_TOP;

% cd data
% save WANG_results_liat mean_corr var_corr f var_f
