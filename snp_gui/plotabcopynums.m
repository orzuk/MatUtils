% Plot two allele copy numbers
% We show the intensity levels of two alleles along with their jointly
% determined genotype calls
function PlotABCopyNums(copy_A, copy_B, genotype_vec, title_str, new_figure_flag)


genotypes = unique(genotype_vec); colorvec = 'bgrkmcxo:.';
[AA, AB, BB, NoCall, call_cell] = genotype_call_into_num();

if(nargin < 5)
    new_figure_flag = 1;
end
if(new_figure_flag)
    figure; 
end


hold on;
% Old - Generic
% for i=1:length(genotypes) % assume there are less then ~8
%     I = find(genotype_vec == genotypes(i));
%     plot(copy_A(I), copy_B(I), ['.' colorvec(i)]);
% end

% New - Genotypes 
legend_vec = {}; ctr=1;
for i= [AA, AB, BB, NoCall]
    I = find(genotype_vec == i);
    if(~isempty(I))
        legend_vec{ctr} = call_cell{i};
        plot(copy_A(I), copy_B(I), ['.' colorvec(i)]);
        ctr=ctr+1;
    end
end
    
xlabel('A copy'); ylabel('B copy');
legend(legend_vec); title(title_str);

