function outStruct = snpquartets(psvals,cdfstruct)
% Create a table of quartets from the output of probesetvalues
%
%   QSTRUCT = SNPQUARTETS(PSVALS,CDFSTRUCT) creates a struct of SNP results
%   where PSVALS is the output of probesetvalues and CDFSTRUCT is the
%   corresponding CDF file structure.
%
%   Example:
%
%   CELStruct=affyread('74A_Hind.CEL');
%   CDFStruct=affyread('Mapping50K_Hind240.cdf');
%   psvals = probesetvalues(CELStruct,CDFStruct,'SNP_A-1712762');
%   probeStruct = snpquartets(psvals,CDFStruct)
%   probeStruct.Quartet(1)
%
%   See Also AFFYREAD, PROBESETVALUES.

% Copyright 2007 The MathWorks, Inc.
% $Revision: 1.1 $   $Date: 2007/05/08 17:47:34 $

numProbePairs = size(psvals,1);
usedMask = false(numProbePairs,1);
PMCol = 7;
MMCol = 14;
groupCol = 19;
dirCol = 20;
numGroups = numel(unique(psvals(:,groupCol)));

quartetCount = 1;
outStruct.ProbeSet = cdfstruct.ProbeSets(psvals(1)+1).Name;
% loop over all entries
for count = 1:numProbePairs
    % figure out if this value has already been used
    if usedMask(count)
        continue
    end
    % make room for the quartet
    theQuartet = zeros(numGroups,2);
    groupCount = 1;
    theQuartet(groupCount,1) = psvals(count,groupCol);
    theQuartet(groupCount,2) = count;
    usedMask(count) = true;
    % now find the corresponding pairs
    for inner = count+1:numProbePairs
        if ~usedMask(inner) && isempty(find(psvals(inner,groupCol)==theQuartet(:,1))) %~ismember(psvals(inner,groupCol),theQuartet(:,1))
            groupCount = groupCount+1;
            theQuartet(groupCount,1) = psvals(inner,groupCol);
            theQuartet(groupCount,2) = inner;
            usedMask(inner) = true;
        end
    end
    % set up the outputs first time through
    if quartetCount == 1
        outStruct.AlleleA = cdfstruct.ProbeSets(psvals(1)+1).GroupNames{1};
        outStruct.AlleleB = char(setxor(cdfstruct.ProbeSets(psvals(1)+1).GroupNames,outStruct.AlleleA));
        outStruct.Quartet = struct('A_Sense_PM',[],'B_Sense_PM',[],...
            'A_Sense_MM',[],'B_Sense_MM',[],...
            'A_Antisense_PM',[],'B_Antisense_PM',[],...
            'A_Antisense_MM',[],'B_Antisense_MM',[]);
    end
    % Fill out a table of the quartet
    for outLoop = 1:sum(theQuartet(:,1)~=0)
        if isequal(cdfstruct.ProbeSets(psvals(1)+1).GroupNames{theQuartet(outLoop,1)},outStruct.AlleleA)
            if psvals(theQuartet(outLoop,2),dirCol) == 1
                outStruct.Quartet(quartetCount).A_Sense_PM = psvals(theQuartet(outLoop,2),PMCol);
                outStruct.Quartet(quartetCount).A_Sense_MM = psvals(theQuartet(outLoop,2),MMCol);
            else
                outStruct.Quartet(quartetCount).A_Antisense_PM = psvals(theQuartet(outLoop,2),PMCol);
                outStruct.Quartet(quartetCount).A_Antisense_MM = psvals(theQuartet(outLoop,2),MMCol);
            end
        else
            if psvals(theQuartet(outLoop,2),dirCol) == 1
                outStruct.Quartet(quartetCount).B_Sense_PM = psvals(theQuartet(outLoop,2),PMCol);
                outStruct.Quartet(quartetCount).B_Sense_MM = psvals(theQuartet(outLoop,2),MMCol);
            else

                outStruct.Quartet(quartetCount).B_Antisense_PM = psvals(theQuartet(outLoop,2),PMCol);
                outStruct.Quartet(quartetCount).B_Antisense_MM = psvals(theQuartet(outLoop,2),MMCol);
            end
        end

    end
    quartetCount = quartetCount +1;
end
