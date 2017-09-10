% Compute overlap f_star vs. n curve
function [mean_corr,var_corr,f,var_f]=new_f_star_vs_n_with_reshape(data,samples)
cd e:\Gene_Corr\Or

alpha=single(0.012);

% data=gene_expression_log2;
% samples=Prognosis;
% N_TOP=267;
N_genes=size(data,1);
N_samples=size(data,2);
N_TOP=floor(alpha*N_genes);
idx_poor=int16(find(samples==1));
idx_good=int16(find(samples==0));
N_poor=length(idx_poor);
N_good=N_samples-N_poor;
ratio=N_poor/N_samples;
N_iterations=single(10);
% min_i=ceil(1/ratio);
% max_i=floor(N_samples/2);
res=single(4);
min_i=single(10);
max_i=single(48);
group_sizes=[min_i:res:max_i];
num_of_sizes=length(group_sizes);
mean_corr=zeros(N_genes,num_of_sizes,'single');
var_corr=zeros(N_genes,num_of_sizes,'single');
f=zeros(1,num_of_sizes,'single');
var_f=zeros(1,num_of_sizes,'single');
segments=4;
for i=1:num_of_sizes
    group_size=group_sizes(i);
    clear curr_samples;
    N_groups=floor(N_samples/group_size);%number of groups
    N_couples=N_groups*(N_groups-1)/single(2);
    curr_N_poor=floor(group_size*ratio);%number of poor samples in each group
    curr_N_good=floor(group_size*(1-ratio));
    clear S
    hashlamot=group_size-(curr_N_poor+curr_N_good);
    end_poor=N_groups*curr_N_poor;
    end_good=N_groups*curr_N_good;
    end_hashlamot=hashlamot*N_groups;
    ff=[];
    var_ff=[];
    group_num=int16(repmat(1:N_groups,group_size,1));
    group_num=int16(repmat(group_num(:)',N_genes,1));%dtermine the location within the column
    idx_in_group=int16(repmat(1:group_size,N_genes,N_groups));%determine the column
    reshuffled_indices=int32(idx_in_group-1)*int32(N_genes*N_groups)+...
        int32(group_num-1)*int32(N_genes)+repmat((1:int32(N_genes))',1,group_size*N_groups);
    clear   group_num  idx_in_group  corr list
%     save information N_TOP N_couples N_genes N_good N_groups N_iterations N_poor N_samples alpha curr_N_good curr_N_poor end_good end_hashlamot end_poor f ff group_size group_sizes hashlamot i idx_good idx_poor max_i mean_corr min_i num_of_sizes ratio res samples var_ff var_f var_corr
%     clear N_TOP N_couples N_genes N_good N_groups N_iterations N_poor N_samples alpha curr_N_good curr_N_poor end_good end_hashlamot end_poor f ff group_size group_sizes hashlamot i idx_good idx_poor max_i mean_corr min_i num_of_sizes ratio res samples var_ff var_f var_corr
    [val,ind_re]=sort(reshuffled_indices(1:end/4));
    sof=length(ind_re);
    for ii=2:segments
         [val,ind_re((ii-1)*sof+1:ii*sof)]=sort(reshuffled_indices((ii-1)*sof+1:ii*sof));
         ind_re((ii-1)*sof+1:ii*sof)=(ii-1)*sof+ind_re((ii-1)*sof+1:ii*sof);
    end
    ind_re=int32(ind_re);
    clear val
%     corr=zeros(N_genes,N_groups,'single');
%     list=zeros(N_TOP,N_groups,'single');
    for j=1:N_iterations
        j
        x=int16(idx_poor(randperm(N_poor)));
        y=int16(idx_good(randperm(N_good)));
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
        clear x y idx_all
       all_data=data(:,S(:));
        all_data=reshape(all_data(reshuffled_indices(:)),N_genes*N_groups,group_size);
%         [norm_data,mean_genes,norm_genes]=dna_N_new(all_data);
       all_data=dna_N(all_data);
       all_data=reshape(all_data(ind_re),N_genes,N_groups*group_size);
        
        all_samples=reshape(samples(S(:)),N_groups,group_size);
        all_samples=dna_N(all_samples);
        for k=1:N_groups
            corr(:,k)=all_data(:,S(:,k))*all_samples(k,:)';
            [val,ind]=sort(-abs(corr(:,k)));
            clear val
            list(:,k)=int16(ind(1:N_TOP));
            for kk=1:k-1
                ff(end+1)=length(intersect(list(:,k),list(:,kk)));               
            end
        end
%         clear all_data all_samples
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
% cd data
% save WANG_results_liat mean_corr var_corr f var_f