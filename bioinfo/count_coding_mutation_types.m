% Count number of missense, synonymous and stop altering mutations in the
% gentic code 
function [n_synonymous, n_missense, n_stop, MutTable, MutType] = count_coding_mutation_types()

n_synonymous = 0; n_missense = 0; n_stop = 0; 

nt_vec = 'ACGT'; 

MutTable = cell(64, 3, 4); MutType = zeros(64, 3, 4); 
codons = fields(geneticcode()); codons = codons(2:end-1); 
for i=1:64 
    AA = nt2aa(codons{i}); % get AA 
    for j=1:3 % which nucleotide to change
        for k=1:4 % which mutation to apply 
            alt_codon = codons{i}; 
            alt_codon(j) = nt_vec(k);             
            if(~strcmp(codons{i}, alt_codon))
               if(nt2aa(alt_codon) == AA)
                   n_synonymous = n_synonymous+1; MutType(i,j,k) = 0; 
               else
                  if((AA == '*') || (nt2aa(alt_codon) == '*'))
                      n_stop = n_stop + 1;  MutType(i,j,k) = 2; 
                  else
                      n_missense = n_missense + 1;   MutType(i,j,k) = 1; 
                  end
               end
            end
            
            MutTable{i}{j}{k} = [codons{i} '->' alt_codon];
        end
    end
end
    
