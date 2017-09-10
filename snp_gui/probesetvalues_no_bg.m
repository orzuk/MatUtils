function [table,lookup] = probesetvalues_no_bg(celStruct,cdfStruct,ID,varargin)
% PROBESETVALUES extracts probe set values from probe results
%
%   PSVALUES = PROBESETVALUES(CELSTRUCT,CDFSTRUCT,PS) creates a table of
%   values for probe set PS from the probe data in a CEL file structure
%   CELSTRUCT, where PS is a probe set index or probe set name from the CDF
%   library file structure CDFSTRUCT.
%
%   PSVALUES is a matrix with 20 columns and one row for each probe pair in
%   the probe set. The columns correspond to these fields:
%
%     'ProbeSetNumber'
%     'ProbePairNumber'
%     'UseProbePair'
%     'Background'
%     'PMPosX'
%     'PMPosY'
%     'PMIntensity'
%     'PMStdDev'
%     'PMPixels'
%     'PMOutlier'
%     'PMMasked'
%     'MMPosX'
%     'MMPosY'
%     'MMIntensity'
%     'MMStdDev'
%     'MMPixels'
%     'MMOutlier'
%     'MMMasked'
%     'GroupNumber'
%     'Direction'
%
%   The 'UseProbePair' column is for backwards compatibility and is not
%   currently used.
%
%   Example:
%       celStruct = affyread('Ecoli-antisense-121502.cel');
%       cdfStruct = affyread(...
%                    'C:\Affymetrix\LibFiles\Ecoli_ASv2\Ecoli_ASv2.CDF');
%       % get the values for probe set 'argG_b3172_at'
%       psvals = probesetvalues(celStruct,cdfStruct,'argG_b3172_at')
%
%   See also AFFYREAD, CELINTENSITYREAD, PROBELIBRARYINFO, PROBESETLINK,
%   PROBESETLOOKUP, PROBESETPLOT, RMABACKADJ.

%   Affymetrix and NetAffx are registered trademarks of Affymetrix, Inc.

% Copyright 2003-2007 The MathWorks, Inc.
% $Revision: 1.1.12.1.4.9 $   $Date: 2007/04/25 16:38:11 $

if nargin <3
    error('Bioinfo:probesetvalues:NotEnoughInputs',...
        'Not enough input arguments.');
end

% Now check that the CEL and CDF struct are cel and cdf structs

if ~isstruct(celStruct)
    error('Bioinfo:probesetvalues:CelStructNotStruct',...
        'The first input must be a structure.');
end

if  ~isstruct(cdfStruct)
    error('Bioinfo:probesetvalues:CdfStructNotStruct',...
        'The second input must be a structure.');
end

if ~isfield(celStruct,'Name') || ~isfield(celStruct,'ChipType') ||...
        ~isfield(celStruct,'Probes') || isempty(regexpi(celStruct.Name,'.cel$'))
    error('Bioinfo:probesetvalues:BadCelStruct',...
        'The first input must be a structure created by AFFYREAD from a CEL file.');
end
if ~isfield(cdfStruct,'Name') || ~isfield(cdfStruct,'ChipType') ||...
        ~isfield(cdfStruct,'ProbeSets') || isempty(regexpi(cdfStruct.Name,'.CDF$'))
    error('Bioinfo:probesetvalues:BadCdfStruct',...
        'The second input must be a structure created by AFFYREAD from a CDF file.');
end

% Check that the ChipType match
if strcmpi(celStruct.ChipType, cdfStruct.ChipType) == 0
    error('Bioinfo:probesetvalues:ChipTypeMismatch',...
        'The CDF ChipType ''%s'' is not the same as the CEL ChipType ''%s''.',...
        cdfStruct.ChipType,celStruct.ChipType);
end

backFlag = false;

% deal with the various inputs
if nargin > 3
    if rem(nargin,2) == 0
        error('Bioinfo:probesetvalues:IncorrectNumberOfArguments',...
            'Incorrect number of arguments to %s.',mfilename);
    end
    okargs = {'background'};
    for j=1:2:nargin-3
        pname = varargin{j};
        pval = varargin{j+1};
        k = find(strncmpi(pname,okargs,numel(pname)));
        if isempty(k)
            error('Bioinfo:probesetvalues:UnknownParameterName',...
                'Unknown parameter name: %s.',pname);
        elseif length(k)>1
            error('Bioinfo:probesetvalues:AmbiguousParameterName',...
                'Ambiguous parameter name: %s.',pname);
        else
            switch(k)
                case 1  % noback- Undocumented, only used by PROBESETPLOT (DEFAULT=TRUE)
                    backFlag = opttf(pval);
                    if isempty(backFlag)
                        error('Bioinfo:probesetvalues:InputOptionNotLogical',...
                            '%s must be a logical value, true or false.',...
                            upper(char(okargs(k))));
                    end
            end
        end
    end
end

% get the ID
% if ischar(theID)
%     ID = find(strncmp(theID,{cdfStruct.ProbeSets.Name},numel(theID)));
%     if isempty(ID)
%         error('Bioinfo:probesetvalues:UnknownProbeName',...
%             'Unknown probe set name: %s.',theID);
%     elseif length(ID)>1
%         warning('Bioinfo:probesetvalues:AmbiguousParameterName',...
%             'Ambiguous probe set name: %s.',theID);
%         ID = ID(1);
%     end
% else
%     ID = theID;
% end

numCols = cdfStruct.Cols;

% Check that the ID is valid
if ~isscalar(ID)
    error('Bioinfo:probesetvalues:PSIDNotScalar',...
        'The probe set ID must be a scalar value.');
end
if ID > numel(cdfStruct.ProbeSets) || ID < 1 ||  floor(ID)~= ID
    error('Bioinfo:probesetvalues:BadPSID',...
        'The probe set ID (%d) is invalid.',ID);
end
NumPairs = cdfStruct.ProbeSets(ID).NumPairs;
table = zeros(NumPairs,18);
lookup = zeros(NumPairs,2);
thePairs = cdfStruct.ProbeSets(ID).ProbePairs;

PMXCol = strcmp('PMPosX',cdfStruct.ProbeSetColumnNames);
PMYCol = strcmp('PMPosY',cdfStruct.ProbeSetColumnNames);
MMXCol = strcmp('MMPosX',cdfStruct.ProbeSetColumnNames);
MMYCol = strcmp('MMPosY',cdfStruct.ProbeSetColumnNames);
GroupCol = strcmp('GroupNumber',cdfStruct.ProbeSetColumnNames);
DirectionCol = strcmp('Direction',cdfStruct.ProbeSetColumnNames);

IntensityCol = strcmp('Intensity',celStruct.ProbeColumnNames);
StdDevCol = strcmp('StdDev',celStruct.ProbeColumnNames);
PixelsCol = strcmp('Pixels',celStruct.ProbeColumnNames);
OutlierCol = strcmp('Outlier',celStruct.ProbeColumnNames);
MaskedCol = strcmp('Masked',celStruct.ProbeColumnNames);

for inner = 1:NumPairs
    PMX = thePairs(inner,PMXCol);
    PMY = thePairs(inner,PMYCol);
    PMRow = PMY*numCols + PMX +1;
    table(inner,1) = ID-1;
    table(inner,2) = inner-1;
    table(inner,5) = PMX;
    table(inner,6) = PMY;
    lookup(inner,1) = PMRow;
    table(inner,7) = celStruct.Probes(PMRow,IntensityCol);
    table(inner,8) = celStruct.Probes(PMRow,StdDevCol);
    table(inner,9) = celStruct.Probes(PMRow,PixelsCol);
    table(inner,10) = celStruct.Probes(PMRow,OutlierCol);
    table(inner,11) = celStruct.Probes(PMRow,MaskedCol);
    MMX = thePairs(inner,MMXCol);
    MMY = thePairs(inner,MMYCol);
    MMRow = MMY*numCols + MMX + 1;
    lookup(inner,2) = MMRow;
    table(inner,12) = MMX;
    table(inner,13) = MMY;
    table(inner,14) = celStruct.Probes(MMRow,IntensityCol);
    table(inner,15) = celStruct.Probes(MMRow,StdDevCol);
    table(inner,16) = celStruct.Probes(MMRow,PixelsCol);
    table(inner,17) = celStruct.Probes(MMRow,OutlierCol);
    table(inner,18) = celStruct.Probes(MMRow,MaskedCol);
    % backgrounds -- fill in later if necessary
    table(inner,4) = 0;
    table(inner,19) = thePairs(inner,GroupCol);
    table(inner,20) = thePairs(inner,DirectionCol);
end
if backFlag
    [adjVals,zones,backgrounds] = zonebackadj(celStruct,'cdf',cdfStruct,'bgindices',lookup(:,1));
    table(:,4) = backgrounds;
end

