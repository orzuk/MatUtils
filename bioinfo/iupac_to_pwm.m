% Convert a iupac code sequence to a pwm 
% We currently don't support numeric (?)
% 
% Input: 
% iupac_seq - sequence in iupac format 
% 
% Output: 
% pwm - 4xL probabilities matrix
% 
function pwm = iupac_to_pwm(iupac_seq)

if(isnumeric(iupac_seq))
    iupac_seq = int2nt(iupac_seq);
end
iupac_seq = upper(iupac_seq); % make sure everything is upper-case 
L = length(iupac_seq);
pwm = zeros(4,L);

% Start with consensus
pwm(1,iupac_seq == 'A') = 1;
pwm(2,iupac_seq == 'C') = 1;
pwm(3,iupac_seq == 'G') = 1;
pwm(4,iupac_seq == 'T') = 1;


% Now go to couples
pwm([1 2],iupac_seq == 'M') = 0.5; % A C 
pwm([1 3],iupac_seq == 'R') = 0.5; % A G (purine)
pwm([1 4],iupac_seq == 'W') = 0.5; % A T (weak)
pwm([2 3],iupac_seq == 'S') = 0.5; % C G (strong)
pwm([2 4],iupac_seq == 'Y') = 0.5; % C T (pyrimidine)
pwm([3 4],iupac_seq == 'K') = 0.5; % G T 


% Do all triplets
pwm([2 3 4],iupac_seq == 'B') = 0.333333333; % C G T (not A)
pwm([1 3 4],iupac_seq == 'D') = 0.333333333; % A G T (not C)
pwm([1 2 4],iupac_seq == 'H') = 0.333333333; % A C T (not G)
pwm([1 2 3],iupac_seq == 'V') = 0.333333333; % A C G (not T)


% Finally the non-informative 
pwm(:,iupac_seq == 'N') = 0.25;
