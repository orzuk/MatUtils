function [mean_corr,var_corr,f,var_f]=f_star_vs_n_with_reshape(data,samples)

ttt = cputime;

cd E:\Gene_Corr\Or

alpha=0.012;

% data=gene_expression_log2;
% samples=Prognosis;
% N_TOP=267;
N_genes=size(data,1);
N_samples=size(data,2);
N_TOP=floor(alpha*N_genes);
idx_poor=find(samples==1);
idx_good=find(samples==0);
N_poor=length(idx_poor);
N_good=N_samples-N_poor;
ratio=N_poor/N_samples;
N_iterations=10;
% min_i=ceil(1/ratio);
% max_i=floor(N_samples/2);
res=4;
min_i=10;
max_i=48;
group_sizes=[min_i:4:max_i];
num_of_sizes=length(group_sizes);
mean_corr=zeros(N_genes,num_of_sizes);
var_corr=zeros(N_genes,num_of_sizes);
f=zeros(1,num_of_sizes);
var_f=zeros(1,num_of_sizes);
for i=1:num_of_sizes
    group_size=group_sizes(i);
    clear corr curr_samples;
    N_groups=floor(N_samples/group_size);%number of groups
    N_couples=N_groups*(N_groups-1)/2;
    curr_N_poor=floor(group_size*ratio);%number of poor samples in each group
    curr_N_good=floor(group_size*(1-ratio));
    clear S
    hashlamot=group_size-(curr_N_poor+curr_N_good);
    end_poor=N_groups*curr_N_poor;
    end_good=N_groups*curr_N_good;
    end_hashlamot=hashlamot*N_groups;
    ff=[];
    var_ff=[];
    corr=zeros(N_genes,N_groups);
    list=zeros(N_TOP,N_groups);
    for j=1:N_iterations
        j
        x=idx_poor(randperm(N_poor));
        y=idx_good(randperm(N_good));
        if(hashlamot>0)
            idx_all=[x(end_poor+1:end);...
                y(end_good+1:end)];
        S=[reshape(x(1:end_poor),curr_N_poor,N_groups);...
            reshape(y(1:end_good),curr_N_good,N_groups);...
            reshape(idx_all(1:end_hashlamot),hashlamot,N_groups)] ;
        else
        S=[reshape(x(1:end_poor),curr_N_poor,N_groups);...
            reshape(y(1:end_good),curr_N_good,N_groups)];
        end
        
        for k=1:N_groups
            norm_data=dna_N(data(:,S(:,k)));
            curr_samples=samples(S(:,k));
            norm_samples=curr_samples-mean(curr_samples);
            norm_samples=norm_samples/sqrt(sum(norm_samples.^2));
            corr(:,k)=norm_data*norm_samples;
            [val,ind]=sort(-abs(corr(:,k)));
            list(:,k)=ind(1:N_TOP);
            for kk=1:k-1
                ff(end+1)=length(intersect(list(:,k),list(:,kk)));               
            end
        end
        mean_corr(:,i)=mean_corr(:,i)+mean(corr,2);
        var_corr(:,i)=var_corr(:,i)+mean(corr.^2,2)-mean(corr,2).^2;
    end
    f(i)=mean(ff);
    var_f(i)=mean(ff.^2)-f(i)^2;
end
mean_corr=mean_corr/N_iterations;
var_corr=var_corr/N_iterations;
f=f/N_TOP;    
var_f=sqrt(var_f)/N_TOP;

cputime - ttt

% cd data
% save WANG_results_liat mean_corr var_corr f var_f