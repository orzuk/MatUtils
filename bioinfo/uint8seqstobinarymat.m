% Transfer the sequence into a binary matrix form. 
% This is a most wastefull format, but enables one to perform some linear
% algebra with the sequences. 
% Input seqs in uint8 form where each base is represented by one
% uint8 number.
% Output a binary matrix 

function BinSeqs = Uint8SeqsToBinaryMat(Seqs)

NumSeqs = size(Seqs, 1); SeqsLen = size(Seqs, 2);



for i=1:4 % loop over sequences
    I{i} = find(Seqs' == i)-1; 
    J{i} = floor(I{i} ./ SeqsLen)+1;
    I{i} = mod(I{i}, SeqsLen);
    I{i} = I{i} * 4 + i;
end


I = [I{1}' I{2}' I{3}' I{4}']; 
J = [J{1}' J{2}' J{3}' J{4}'];
%% BinSeqs  = (single(full(sparse(I, J, ones(1, length(I)),  4*SeqsLen, NumSeqs))))';
BinSeqs  = (double(accumarray([I', J'], ones(1, length(I))',  [4*SeqsLen, NumSeqs])))';
