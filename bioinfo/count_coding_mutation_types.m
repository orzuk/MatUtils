% Count number of missense, synonymous and stop altering mutations in the gentic code 
% Output: 
% n_synonymous - number of synonymous mutations 
% n_synonymous - number of synonymous mutations 
% n_synonymous - number of synonymous mutations 
% MutTable - Table with all mutations, with codons and amino acid changes
% MutType - Table with types of all mutations (0 - synonymous, 1 - missense, 2 - stop)
% AAMutTable - number of mutations for each type for every amino acid
%
function [n_synonymous, n_missense, n_stop, MutTable, MutType, AAMutTable] = count_coding_mutation_types()

n_synonymous = 0; n_missense = 0; n_stop = 0; 
nt_vec = 'ACGT'; 

MutTable = cell(64, 3, 4); MutType = zeros(64, 3, 4); 
codons = fields(geneticcode()); codons = codons(2:end-1); 
AAMutTable = zeros(25, 3); 
for i=1:64 % loop on codons 
    AA = nt2aa(codons{i}, 'AlternativeStartCodons', false); % get AA 
    for j=1:3 % which nucleotide to change
        for k=1:4 % which mutation to apply 
            alt_codon = codons{i}; 
            alt_codon(j) = nt_vec(k);             
            MutTable{i}{j}{k} = [codons{i} '->' alt_codon ', ' AA '->' nt2aa(alt_codon, 'AlternativeStartCodons', false)];
            if(~strcmp(codons{i}, alt_codon)) % get rid of 'same' mutation 
               if(nt2aa(alt_codon, 'AlternativeStartCodons', false) == AA) % synonymous
                   n_synonymous = n_synonymous+1; MutType(i,j,k) = 0; 
               else 
                  if((AA == '*') || (nt2aa(alt_codon, 'AlternativeStartCodons', false) == '*')) % stop 
                      n_stop = n_stop + 1;  MutType(i,j,k) = 2; 
                  else
                      n_missense = n_missense + 1;   MutType(i,j,k) = 1; 
                  end
               end
               AAMutTable(aa2int(AA), MutType(i,j,k)+1) = AAMutTable(aa2int(AA), MutType(i,j,k)+1)+1;
            end
        end
    end
end
AAMutTable = num2cell(AAMutTable); 
for i=1:25
    AAMutTable{i,4} = int2aa(i);     AAMutTable{i,5} = aminolookup(int2aa(i)); 

end
AAMutTable = [{'Synonymous', 'missense', 'stop', 'AA', 'AA-Full'}' AAMutTable']'; 

