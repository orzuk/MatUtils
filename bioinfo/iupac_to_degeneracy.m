% Convert a iupac code sequence to degeneracy
% 
% Input: 
% iupac_seq - sequence in iupac format 
% 
% Output: 
% degeneracy - how many different bases are in each position
% 
function degeneracy = iupac_to_degeneracy(iupac_seq)

if(isnumeric(iupac_seq))
    iupac_seq = int2nt(iupac_seq);
end
iupac_seq = upper(iupac_seq); % make sure everything is upper-case 
L = length(iupac_seq);
degeneracy = zeros(1,L);

degeneracy(regexp(iupac_seq, '[ACGT]')) = 1;  % Start with consensus
degeneracy(regexp(iupac_seq, '[MRWSYK]')) = 2; % Now go to couples
degeneracy(regexp(iupac_seq, '[BDHV]')) = 3; % Do all triplets
degeneracy(regexp(iupac_seq, '[N]')) = 4; % Finally the non-informative 
