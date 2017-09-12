%function [AA, AB, BB, NoCall, call_cell] = genotype_call_into_num()
function [AA, AB, BB, NoCall, call_cell] = genotype_call_into_num()

AssignAllGlobalConstants();  %% AA = 1; AB = 2; BB = 3; NoCall = 4;

call_cell = {'AA';'AB';'BB';'NoCall'};