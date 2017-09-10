% Collapse a genotype count table for case/control study according to test
% at hand
%
% Input:
% allele_tab - a 2X2 table of allelic counts
% test_str - string saying which test we're using
% data_format - whether to use 2X2 and 2X3 tables (default) or 4XN and 6XN
%
% Output:
% genotype_tab - a 2X3 table of genotype counts
%
function genotype_tab = expand_allele_table(allele_tab, test_str, data_format)

if(~exist('data_format', 'var') || isempty(data_format))
    data_format = '2X2';
end
switch data_format
    case '2X2'
        genotype_tab = zeros(3,2);
    case '4XN'
        genotype_tab = zeros(size(allele_tab,1), 6);
end


if(~exist('test_str', 'var') || isempty(test_str))
    test_str = 'armitage';
end
switch test_str
    case {'trend', 'additive', 'armitage'} % assume HW equilibrium
        switch data_format
            case '2X2'
                p_case_controls = [sum(allele_tab(:,1)) sum(allele_tab(:,2))];  % probability of cases and of controls
                for i=1:2
                    genotype_tab(:,i) = [allele_tab(1,i)^2 ...
                        2*allele_tab(1,i)*allele_tab(2,i) allele_tab(2,i)^2] ./ p_case_controls(i);
                end
            case '4XN'
                p_case_controls = [sum(allele_tab(:,[1 3]),2) sum(allele_tab(:,[2 4]),2)];  % probability of cases and of controls
                for i=1:2 % loop on disease status
                    genotype_tab(:,([1 3 5]) + (i-1)) = [allele_tab(:,i).^2 ...
                        2.*allele_tab(:,i).*allele_tab(:,i+2) allele_tab(:,i+2).^2] ./ ...
                        repmat(p_case_controls(:,i), 1, 3);
                end
        end
    case {'dominance', 'dominant'}
        
    case 'recessive'
        
end

