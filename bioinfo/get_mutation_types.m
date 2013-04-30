% Compute a 64x64 table with all possible mutation types
% 
% Output: 
% mutation_types_table - a 64x64 table saying for each substitution its type
% codons - list of all 64 codons  
%
function [mutation_types_table codons] = get_mutation_types()

Assign24MammalsGlobalConstants;

mutation_types_table = zeros(64); 
codons = int2nt(unpack_seqs((0:63)', 3)); codons = sortrows(codons);
amino_acids = nt2aa(codons); 

for i=1:64    
    for j=1:64
        if(amino_acids(i) == amino_acids(j))
            mutation_types_table(i,j) = SYNONYMOUS;
        else
            mutation_types_table(i,j) = MISSENSE;
            if(amino_acids(i) == '*') % stop codon lost
                mutation_types_table(i,j) = STOP_LOST;
            end
            if(amino_acids(j) == '*') % stop codon gained
                mutation_types_table(i,j) = STOP_GAINED;
            end                
        end
    end
end




