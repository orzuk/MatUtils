% Written by Or Zuk 6/2007
%
% The function saves the data that we need to display after the HMM,
% in the display format to be read by the chromosome viewer.
%
% Output the structure DispStruct, with the following fields:
% SampleName - The name of the sample
% ChipName - The name of the chip
% GenomeBuild - The genome build version used
% Chrom - One for each chromosome:
% Chrom{i}.Locs - The SNPs locations on chromosome #i
% Chrom{i}.SNPsIDs - The SNPs IDs on chromosome #i
% Chrom{i}.Data - The SNPs intensity/data values/whatever for chromosome
%                 #i. We use Three data points for each SNP - the alpha copy
%                 and beta copy, which are discrete, and then the total copy which is continuous !!!  (but one can use any number)
% Chrom{i}.Segments - The segments we want to highlight. There are 3 numbers
% for each SNP. The format is: [start, end, value].
% Note that we do not want to create segments which go beyond p-q arm (centromerss)
% Chrom{i}.Genotypes - This is new - we keep here the genotypes:
% 0 - AA, 1-AB, 2-BA 3-BB. Currently 0 and 3 should be displayed the same
% and 1 and 2 also. (only tell hetro/homo), but this may change in the future.
function  [display_file_name DispStruct] = DisplayFormatSave(user_dir, sample_name, chip_type, genome_assembly, ...
    CHROM_MATS, HMMParamsStruct, VitStruct, ProbsStruct)

Upperbound = 2.5; % above this we say it's amplifications
Lowerbound = 1.5; % below this we say it's deletions
display_dir = [user_dir  '\display'];
display_file_name = fullfile('display',[sample_name '_' chip_type '_disp']);
DispStruct = [];
DispStruct.SampleName = sample_name;
DispStruct.ChipName = chip_type;
DispStruct.GenomeBuild = genome_assembly;
end_p_location = load_end_p_loc(); % Get the splitting between P and Q
for i=1:HMMParamsStruct.ChromosomesToRun
    DispStruct.Chrom{i} = [];
end
for i=HMMParamsStruct.ChromosomesToRun
    if(~isempty(VitStruct{i}))
        DispStruct.Chrom{i}.Data = [VitStruct{i}.alpha_copy VitStruct{i}.beta_copy];
        % calc total raw copy num   (fractional value !)
        DispStruct.Chrom{i}.Data(:,3) = sum(ProbsStruct{i}.total_copy .* repmat([0:4], size(ProbsStruct{i}.total_copy, 1), 1), 2);
        % calc RAW total raw copy num   (fractional value !)
        DispStruct.Chrom{i}.Data(:,4) = smooth( CHROM_MATS.data_A{i} + CHROM_MATS.data_B{i}, HMMParamsStruct.SmoothWindow ); % raw intensity (after normalization)
        DispStruct.Chrom{i}.Genotypes = VitStruct{i}.joint_genotype;
        DispStruct.Chrom{i}.Locs = CHROM_MATS.loc{i};
        DispStruct.Chrom{i}.SNPsIDs = CHROM_MATS.snp_ids{i}; % added also SNPs IDs

        AmpInds = (VitStruct{i}.alpha_copy+VitStruct{i}.beta_copy > Upperbound)'; % Find the segments positions
        [AmpIndsStarts AmpIndsEnds] = GetSegmentsFromIndexes(AmpInds);
        DelInds = (VitStruct{i}.alpha_copy+VitStruct{i}.beta_copy < Lowerbound)';
        [DelIndsStarts DelIndsEnds] = GetSegmentsFromIndexes(DelInds);
        AllIndsStarts = [AmpIndsStarts DelIndsStarts];
        AllIndsEnds = [AmpIndsEnds DelIndsEnds];
        [AllIndsStarts SortPerm] = sort(AllIndsStarts); AllIndsEnds = AllIndsEnds(SortPerm); % Sort segments


        if(~isempty(AllIndsStarts))
            % Now we must split the segments which fall on both P and Q arms:
            SegmentsToSplit = intersect( find(CHROM_MATS.loc{i}(AllIndsStarts) < end_p_location(i)), ...
                find(CHROM_MATS.loc{i}(AllIndsEnds) > end_p_location(i)) );
            if(~isempty(SegmentsToSplit))
                last_p_ind = max(find(CHROM_MATS.loc{i} < end_p_location(i)));
                first_q_ind = min(find(CHROM_MATS.loc{i} > end_p_location(i)));
                for j=1:length(SegmentsToSplit)
                    AllIndsStarts = [AllIndsStarts(1:SegmentsToSplit(j)) first_q_ind  AllIndsStarts(SegmentsToSplit(j)+1:end)];
                    AllIndsEnds = [AllIndsEnds(1:SegmentsToSplit(j)-1)  last_p_ind AllIndsEnds(SegmentsToSplit(j):end)];
                end
            end

            % Now get rid of empty SNPs
            EmptyIndsStarts = find(strcmp('', CHROM_MATS.snp_ids{i}(AllIndsStarts)));
            EmptyIndsEnds = find(strcmp('', CHROM_MATS.snp_ids{i}(AllIndsEnds)));

            for j=1:length(EmptyIndsStarts)
                while(   (isempty(CHROM_MATS.snp_ids{i}{AllIndsStarts(EmptyIndsStarts(j))})) & ...
                        (AllIndsStarts(EmptyIndsStarts(j)) <= length(CHROM_MATS.snp_ids{i})) )
                    AllIndsStarts(EmptyIndsStarts(j)) = AllIndsStarts(EmptyIndsStarts(j))+1;
                end
            end
            for j=1:length(EmptyIndsEnds)
                while(   (isempty(CHROM_MATS.snp_ids{i}{AllIndsEnds(EmptyIndsEnds(j))})) & ...
                        (AllIndsEnds(EmptyIndsEnds(j)) >= 1)  )
                    AllIndsEnds(EmptyIndsEnds(j)) = AllIndsEnds(EmptyIndsEnds(j))-1;
                end
            end

            % Remove too short and empty SNPs segments
            BadSegments = find(AllIndsStarts == AllIndsEnds);
            BadSegments = union(BadSegments, find(strcmp('', CHROM_MATS.snp_ids{i}(AllIndsEnds))));
            BadSegments = union(BadSegments, find(strcmp('', CHROM_MATS.snp_ids{i}(AllIndsStarts))));
            GoodSegments = setdiff(1:length(AllIndsStarts), BadSegments);
            AllIndsStarts = AllIndsStarts(GoodSegments);
            AllIndsEnds = AllIndsEnds(GoodSegments);

            if(~isempty(AllIndsStarts))
                DispStruct.Chrom{i}.Segments(:,1) = CHROM_MATS.loc{i}(AllIndsStarts);
                DispStruct.Chrom{i}.Segments(:,2) = CHROM_MATS.loc{i}(AllIndsEnds);
                V_VEC_cumsum = [0 cumsum(VitStruct{i}.alpha_copy+VitStruct{i}.beta_copy)'];
                DispStruct.Chrom{i}.Segments(:,3) = (V_VEC_cumsum(AllIndsEnds+1) - V_VEC_cumsum(AllIndsStarts)) ./ ...
                    (AllIndsEnds-AllIndsStarts+1); % Find the average copy # of each semgment
            else
                DispStruct.Chrom{i}.Segments = [];  % Put an empty segment
            end
        else
            DispStruct.Chrom{i}.Segments = [];  % Put an empty segment
        end % If notempty
    else % If notempty
        DispStruct.Chrom{i} = [];  
    end

end
if(~exist(display_dir)) % Make the directory if needed
    eval(['mkdir ' display_dir]);
end
save(fullfile(user_dir, display_file_name),'DispStruct'); % Save the display structure
