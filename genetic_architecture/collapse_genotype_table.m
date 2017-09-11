% Collapse a genotype count table for case/control study according to test
% at hand
% 
% Input: 
% genotype_tab - a 2X3 table of genotype counts 
% test_str - how to collaps (which test)
% 
% Output: 
% allele_tab - a 2X2 table of allele counts 
% 
function allele_tab = collapse_genotype_table(genotype_tab, test_str)

allele_tab = zeros(2); 
switch test_str
case {'trend', 'additive', 'armitage'}
    allele_tab = [2.*genotype_tab(1,:) + genotype_tab(2,:); ...
        genotype_tab(2,:) + 2.*genotype_tab(3,:)] ./ 2; 
    case {'dominance', 'dominant'}
        
        case 'recessive'
            
end

