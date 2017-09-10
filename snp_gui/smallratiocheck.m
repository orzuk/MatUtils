% A small check to see that two .txt files (one with genotypes and one with probe intensities match 
ShortGenotype = loadcellfile('E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\ShortGenotypesOneSample.txt');
ShortProbeCalls = loadcellfile('E:\Research\HMM_Chromosome\SNP\HAPMAP\AffyBenchmark_data\ShortNA06985_Hind_B5_3005533.txt');

ShortGenotypeSnpIDs = ShortGenotype(3:end,2); ShortProbeCallsSnpIDs = ShortProbeCalls(3:end,2);

[snp_intersect I J] = intersect(ShortGenotypeSnpIDs, ShortProbeCallsSnpIDs);

should_be_zero = max( (I-J).^2)

ShortGenotypeNum = zeros(1,length(ShortGenotype(3:end,:)));
ShortGenotypeNum(strmatch('AA', ShortGenotype(3:end,3))) = 0;
ShortGenotypeNum(strmatch('AB', ShortGenotype(3:end,3))) = 1;
ShortGenotypeNum(strmatch('BB', ShortGenotype(3:end,3))) = 3;
ShortGenotypeNum(strmatch('NoCall', ShortGenotype(3:end,3))) = 4;

for j=1:size(ShortProbeCalls,2)
    fff = strmatch('null', ShortProbeCalls(:,j));
    if(~isempty(fff))
        for i=fff'
            ShortProbeCalls{i,j} = '1000';
        end
    end
end

ShortProbeCallsNum = str2num(char(ShortProbeCalls(3:end,3:end)));
ShortProbeCallsNum = reshape(ShortProbeCallsNum, 3000, 56);

PA_inds = zeros(1,56); PB = zeros(1,56);
for j=1:56
    if(~isempty(findstr('PA', ShortProbeCalls{2,j+2})))
        PA_inds(j) = 1;
    end
    if(~isempty(findstr('PB', ShortProbeCalls{2,j+2})))
        PB_inds(j) = 1;
    end

end
PA_inds = find(PA_inds);
PB_inds = find(PB_inds);

    ShortAlleleRatio = mean(ShortProbeCallsNum(:,PA_inds),2) ./  mean(ShortProbeCallsNum(:,PB_inds),2);
ShortCopyNum = 2+0.1*randn(3000,1);
ShortCopyNum = mean(ShortProbeCallsNum(:,PA_inds),2) +  mean(ShortProbeCallsNum(:,PB_inds),2);
PlotAlleleRatios(ShortCopyNum, ShortAlleleRatio, ShortGenotypeNum, 'Short Check');
%LOOKS GOOD !!!!! 08/04/07 !!! 



