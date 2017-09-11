% Compute phantom heritability based on values from two science papers
function [V V_add V_add2 h_phantom] = compute_bacterial_phantom_heritability(paper_str, RAF)

if(~exist('RAF', 'var') || isempty(RAF))
    RAF = 0.5; 
end

switch paper_str    
    case 'Chou'
        GA_tab = log([1 1.166 1.509 1.639 ...
            1.096 1.299 1.614 1.812 ...
            1.142 1.320 1.623 1.784 ...
            1.281 1.435 1.752 1.935])'; 
        N=4;
    case 'Khan'
        N=5;
        R = ReadDataFile('../../common_disease_model/data/bacteria/khan_data.txt');
        GA_tab = log(R.Relative_Fitness); % Take the fitness (need to permute!!!)
        mutations_str = unique(cell2vec(R.Genotype));
        genotype_vec = zeros(2^N,N); 
        for j=1:length(mutations_str) 
           genotype_vec(strfind_cell(R.Genotype, mutations_str(j)),j) = 1; 
        end
end

%coordinate_vec = 0:2^N-1
beta = zeros(N,1); 
for i=1:N% Compute narrow sense heritability
    beta(i) = corr(GA_tab, genotype_vec(:,i)); % bitget(0:2^N-1, i)');   
    tmp_beta = polyfit(genotype_vec(:,i), GA_tab, 1);   
    beta2(i) = tmp_beta(1);  
end
V_add = sum(beta.^2 .* RAF .* (1-RAF)); % get V additive 
V_add2 = sum(beta2.^2 .* RAF .* (1-RAF)); % get V additive 

V = var(GA_tab); % assume RAF=0.5
h_phantom = 1 - V_add2 / V;

