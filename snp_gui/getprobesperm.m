% Compute a permutation on 56 elements. The output should be
% in the following order:
% Perfect Match A (14) | Perfect Match B (14) | Miss Match A (14) | Miss Match A (14)
% Inside each set of 14 the ordering is by sense and then by quartets:
% Sense (Q1, Q2, .., Q7) | AntiSense(Q1,Q2,...,Q7)
function ProbesPerm = GetProbesPerm(probes);



Quartets = []; % which quartet (from 1 to 7)
PM = []; % 1 - Perfect Match 0 - Miss Match
AB = []; % 1 - A 0 - B
Sense = []; % 0 - Sense  1 - Antisense
for i=1:56
    Quartets(i) = str2num(probes{i}(2));
    PM(i) = (probes{i}(end-1) == 'P');
    AB(i) = (probes{i}(4) == 'A');
    Sense(i) = (probes{i}(6) == 'A');
end

PM_Inds = find(PM);  MM_Inds = find(PM == 0);
A_Inds = find(AB); B_Inds = find(AB == 0); 

PMA_Inds = intersect(PM_Inds, A_Inds);
PMB_Inds = intersect(PM_Inds, B_Inds);
MMA_Inds = intersect(MM_Inds, A_Inds);
MMB_Inds = intersect(MM_Inds, B_Inds);


ProbesPerm = [PMA_Inds PMB_Inds MMA_Inds MMB_Inds];

