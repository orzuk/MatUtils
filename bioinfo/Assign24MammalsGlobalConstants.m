% All kinds of constants and global variables for the 24 mammals project
Assign24MammalsBasicGlobalConstants; % assign basic constants. From now on call to make defaults.
[machine machine_delim html_outdir] = get_machine_type();

if(machine == UNIX) 
    fig_format_vec = {'fig'}; 
end

% Note: setting default organism_str and genome_ver onlt if they don't exist.
% If the user wants different ones he MUST verify that they exist in the
% code whenever this script is called!
% organism_str = 'HUMAN'; genome_version = 'hg18'; % version of human genome
if(~exist('genome_version', 'var'))
    genome_version = 'hg18'; %% 'hg17'; % 'mm9'; % mm9 is mouse genome (the latest version, matching the alignment)
end
if(isempty(genome_version)) % also enable to give empty genome versions
    genome_version = 'hg18';
end
if(~exist('beta', 'var')) % The 'temperature' (giving the relative importance of strong vs. weak binding sites)
    beta = -2;
end

% check if organism string doesn't match the genome version and adjust
% according to genome version
% if(~exist('organism_str', 'var'))
organism_str = 'MOUSE'; % here we choose default organism to work with
for i_org_str=1:length(genome_versions_vec)
    if(~isempty(strmatch(genome_version, genome_versions_vec{i_org_str}{2}, 'exact')))
        organism_str = genome_versions_vec{i_org_str}{1};
        break; % break loop
    end
end
% end
stat_flag_vec = { 'pi_lods', 'pi_KL', 'branch_length', 'omega', 'omega_lods', 'parsimony'}; % , 'omega'}; %  'pi_lods', 'pi_KL', 'branch_length'}; % 'pi_lods'; % which statistic to use. Current options are: 'pi_lods', 'pi_KL', 'omega', 'branch_length'
switch machine
    case UNIX
        mammals_data_dir = '/seq/orzuk/24mammals/data';
    case PC
        if(~isempty(strfind(lower(pwd), 'c:')))
            mammals_data_dir = 'c:\research\24mammals\data';
        else
            mammals_data_dir = 'T:\24mammals\data';
        end
        pc_additional_mammals_data_dir = 'Z:\24mammals\data';
end
if(exist('pwms_dir', 'var')) % save previous matrix
    save_pwms_dir = pwms_dir;
end

local_genome_dir = fullfile('../data', genome_version); % temp location of genome sequence file 
switch genome_version % organism_str
    case 'mm9'
        %        pwms_dir = '/ahg/regev/mguttman/24Mammals/Pi'; % Here are the input pi values for mouse from Manuel
        pwms_dir = '/seq/orzuk/24mammals/data/mm9/Pi'; % New! Pi data has moved! (this contains only the .mat data!)
        mammals_alignment_dir = '/seq/rinnscratch/mguttman/24Mammals/multiz30Way/ '; % '/ahg/regev/mguttman/24Mammals/multiz30Way'; % the mouse mm9-centered alignment maf format - Erased! Need to change it!!!
        genome_regions_file = '../data/mm9/mm9Annotations.mat';
    case 'hg17' % old 21-way alignment
        pwms_dir = '/ahg/scr3/genome/2x_aln/estimation/mammals21/runs/pi'; % Here are the input pi values from Manuel
        omega_dir = '/ahg/scr3/genome/2x_aln/estimation/mammals21/runs/sampling_four'; % Here are the input omega values from Manuel
        mammals_alignment_dir = '/ahg/scr3/genome/2x_aln/filtered_alignments'; % the human hg17-centered alignment fasta format
        genome_regions_file = '../data/hg17/hg17Annotations.mat';
    case 'hg18' % new 30-way alignment
        %        pwms_dir = '/ahg/scr3/genome/2x_aln/30mammals/eutherian/pi/'; % The new 30-way alignment (hg18 coordinates)
        pwms_dir = '/seq/mgscratch/30mammals/eutherian/pi';
%        pwms_dir = '/seq/mgscratch/30mammals/eutherian/pi' % here put the        /ahg/ location from Evan
        %        pwms_dir = '/seq/orzuk/24mammals/data/hg18/Pi'; % New! Pi data has moved! (this contains only the .mat data!)
        % omega_dir = '/seq/mgscratch/30mammals/eutherian/omega'; % Temp location
        omega_dir = '/ahg/scr4/30mammals/eutherian/omega/'; % permanent location
        omega_12mer_dir = '/seq/mgscratch/30mammals/eutherian/omega/12mer/'; % location for 12mers
        %  omega_dir = '/ahg/scr3/genome/2x_aln/estimation/mammals21/runs/sampling_four'; % Here are the input omega values from Manuel
        mammals_alignment_dir = '/ahg/scr3/mammals/ucsc/multiz44way'; % the human hg18-centered alignment maf format
        genome_regions_file = '../data/hg18/hg18Annotations.mat'; % Note: We still don't have 18 annotations !!!!
        mammals_motifs_dir = '/ahg/scr4/30mammals/eutherian/motifs/'; % sometimes also pwms
        mammals_pi_dir = '../data/hg18/Pi/'; % my PI values
        mammals_model_dir = '/seq/mgscratch/30mammals/models/'; % this one doesn't work:  '/ahg/scr3/genome/2x_aln/30mammals/models/';
        mammals_model_file = fullfile(mammals_model_dir, 'autosomal.mod'); % autosomal model (there's also chr specific model)
        four_fold_dir = '/seq/mgscratch/30mammals/4Dsites'; % /ahg/scr4/30mammals
        masked_genome_regions_file = ...
            '/seq/orzuk/24mammals/data/hg18/masked_regions/all.mask.human.hg18.mat'; % file with regions which are not alignable

    case 'hg19' % New! enable new genome version 
        mammals_alignment_dir = '/ahg/scr4/hg19/maf';
        phylop_dir = '/ahg/scr4/30mammals/phylop';  % new: add phylop conservation directory 

        
    otherwise
        pwms_dir =  '/ahg/scr3/genome/2x_aln/estimation/mammals21/runs/pi';  % same as human for now
end
if(exist('save_pwms_dir', 'var')) % reload previous matrix
    pwms_dir = save_pwms_dir;
end
genome_dir = '/seq/genome/';
encode_dir = '../data/hg18/ENCODE';
encode_regions_file = fullfile(encode_dir, 'encodeRegions.mat');
encode_pwms_dir = fullfile(encode_dir, 'Pi');
encode_omega_dir = fullfile(encode_dir, 'omega'); 
encode_alignment_dir = '/seq/mgscratch/ENCODE_freeze_01_09/ALN_JAN-2009/'; % encode alignments
pwms_HMRD_dir = '/ahg/scr4/30mammals/HMRD/pi/'; % pi files ready and computed
omega_HMRD_dir = '/ahg/scr4/30mammals/HMRD/omega/'; % no omega files ready!!! need to compute them!!!
local_omega_HMRD_dir = '../data/hg18/omega/sub_clade_41971377'; % here's a local version with omega files 

local_mammals_alignment_dir = fullfile('../data/', genome_version, 'alignment'); % get the local path
local_pwms_dir = fullfile('../data/', genome_version, 'Pi'); % get the local path
local_omega_dir = fullfile('../data/', genome_version, 'omega'); % get the local path
local_four_fold_dir = fullfile('../data', genome_version, '4D'); % get the local path
local_model_dir = '../data/trees'; % my local trees
local_model_file = fullfile(local_model_dir, 'autosomal.mod');

mammals_tree_file = '../docs/trees/mammals44tree.txt'; % The phylogenetic tree we use (need to check if it's correct!!!)
% mammals_model_file = '/ahg/scr4/30mammals/models/chr3_4D.mod'; % phylogenetic model (should be identical to the above in the future!!!)
mammals_model_file = '../data/trees/autosomal_full.txt'; % manuel_eutherian_tree.txt'; %  A local copy of Tim's tree (not all 44!)

mammals_46way_model_file = '/seq/orzuk/24mammals/data/trees/vertebrate46way.mod'; % human-based (hg19) 46way alignment
mammals_60way_model_file = '/seq/genome/mouse/mouse_Mm10/mm10.60way.phastCons.mod'; % mouse-based (mm10) 60way alignment



tarjei_chip_dir = ['/seq/tarjei02/SolexaChIP/Sampled_bestprox2_6MM_25.' genome_version]; % directory with many raw solexa chip-seq files
if(~exist('seqs_dir', 'var'))
    seqs_dir = fullfile('../data/', genome_version);
end

% Old way:
% organism_ind = find(strcmp(organism_str, MANUEL_SPECIES_ORDER)); % get from organism string to index
% organism_chr = NUM_SPECIES_CHROMS(organism_ind);
organism_ind = -1;
for i_org_str=1:length(genome_versions_vec)
    if(strcmp(organism_str, genome_versions_vec{i_org_str}{1}))
        organism_ind = i_org_str;
        break;
    end
end
organism_chr = genome_versions_vec{organism_ind}{3};

% genome_annotations_file = fullfile('..', 'data', 'hg17',
% 'hg17Annotations.mat'); % for now only works for hg17, mm9
switch genome_version
    case {'hg17', 'hg18'}
        genome_annotations_file = fullfile('..', 'data', genome_version, [genome_version 'Annotations.mat']); % make it generic (file missing for most genome versions)
    case {'mm8', 'mm9'}
        genome_annotations_file = fullfile('..', 'data', genome_version, [upper(genome_version) 'RefSeq.annotation.mat']); % make it generic (file missing for most genome versions)
end


