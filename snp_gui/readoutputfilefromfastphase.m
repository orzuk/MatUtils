% The function reads the genotypes from the fastPhase output format, 
% and saves them in the matrix M 
function M = ReadOutputFileFromFastPhase(FastPhaseOutputFile, SampleNames)

% First read like an idiot to a cell-file
R = loadCellFile(FastPhaseOutputFile);

Nsamples = size(R,1)/3-7; 
Ngenes = size(str2num(R{22}),2); % number of SNPs

M = zeros(2*Nsamples, Ngenes); % This is the output matrix
for i=1:Nsamples
    M(2*i-1,:) = str2num(R{19+3*i});
    M(2*i,:) = str2num(R{19+3*i+1});
end    
    














