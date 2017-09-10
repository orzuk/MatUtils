% The function takes the genotypes and saves them in a file which is in the
% fastPhase format, to be read by the fastPhase program
function Dummy = PrepareInputFileForFastPhase(data_mat, FastPhaseInputFile, SampleNames)

% Number of SNPs and samples
Nsamples = size(data_mat.data_genotype_A,1);
Ngenes = size(data_mat.data_genotype_A, 2);
FileName = fopen(FastPhaseInputFile, 'wt');

fprintf(FileName, '%d\n',Nsamples);
fprintf(FileName, '%d\n', Ngenes);

% % Loop on samples
for i=1:Nsamples
    fprintf(FileName, '# id %ld  %s\n', i, SampleNames{i});
    for j=1:Ngenes
        if(data_mat.data_genotype_A(i,j) >= 0)
            fprintf(FileName, '%ld', data_mat.data_genotype_A(i,j)); % print 0 or 1
        else
            fprintf(FileName, '?'); % Unknown genotype
        end
    end
    fprintf(FileName, '\n');
    for j=1:Ngenes
        if(data_mat.data_genotype_B(i,j) >= 0)
            fprintf(FileName, '%ld', data_mat.data_genotype_B(i,j)); % print 0 or 1
        else
            fprintf(FileName, '?'); % Unknown genotype
        end
    end
    fprintf(FileName, '\n');
end


fclose(FileName);








