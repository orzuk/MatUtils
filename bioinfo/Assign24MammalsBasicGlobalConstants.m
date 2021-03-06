% All kinds of constants and global variables for the 24 mammals project.
% Here put ONLY basic constants. Do not decide on anything like genome
% version etc.

AssignGeneralConstants(); % Set general constants
AssignStatsConstants(); % set statistical constants

UNKNOWN=-1; INTERGENIC=1; EXON=2; UTR3=3; INTRON=4; UTR5=5; PROMOTER=6; % TSS=7;
REPEAT = 10; AR = 11; ANCESTRAL_REPEATS = AR; SINE = 12; LINE = 13; % repeats
ENHANCER = 20; % numeric values to different genome regions
MIRNA = 30; PIWIRNA = 31; SIRNA = 32; LINCRNA = 33; % non-coding RNAs
FOUR_FOLD = 40; D4 = FOUR_FOLD; TWO_FOLD = 41; SYNONYMOUS = 42; NONSYNONYMOUS = 43; STOP = 44; STOP_GAINED = STOP; STOP_LOST = 45; 
    FRAMESHIFT = 46; MISSENSE = 47; NONSENSE = 48; SPLICESITE = 49; %  sites in proteins
ENCODE = 50; % encode regions 
GENOMIC=100; 

genomic_types_vec = [UNKNOWN INTERGENIC EXON UTR3 INTRON UTR5 PROMOTER ...
    REPEAT ENHANCER MIRNA PIWIRNA SIRNA LINCRNA]; % numeric values to different genome regions
genome_types_synonyms = cell(1, GENOMIC); % allow alternative spelling
genome_types = {'intergenic', 'exon', 'utr3', 'intron', 'utr5', 'promoter'}; % matching order to python file ReadHg17File.py
genome_types_synonyms{EXON} = {'exon', 'exonic'}; 
genome_types_synonyms{INTRON} = {'intron', 'intronic'}; 
genome_types_synonyms{UTR3} = {'utr3', 'utr_3', '3utr', '3_utr'}; 
genome_types_synonyms{UTR5} = {'utr5', 'utr_5', '5utr', '5_utr'}; 
genome_types{REPEAT} = 'repeat';
genome_types{AR} = 'ancestral_repeat';
genome_types{LINE} = 'line';
genome_types{SINE} = 'sine';
genome_types{ENHANCER} = 'enhancer';
genome_types{MIRNA} = 'mirna';
genome_types{PIWIRNA} = 'piwirna';
genome_types{SIRNA} = 'sirna';
genome_types{LINCRNA} = 'lincrna';
genome_types{FOUR_FOLD} = '4d';
genome_types{TWO_FOLD} = '2d';
genome_types{SYNONYMOUS} = 'synonymous'; genome_types_synonyms{SYNONYMOUS} = {'synonymous', 'coding_synonymous'}; 
genome_types{NONSYNONYMOUS} = 'nonsynonymous';
genome_types{MISSENSE} = 'missense';
genome_types{NONSENSE} = 'nonsense'; genome_types_synonyms{NONSENSE} = {'nonsense', 'coding_notmod3'}; 
genome_types{STOP} = 'stop'; genome_types{STOP_GAINED} = 'stop'; genome_types_synonyms{STOP} = {'stop', 'stop_gained'};   % overload stop-gained and stop 
genome_types{STOP_LOST} = 'stop_lost'; 
genome_types{SPLICESITE} = 'splicesite'; genome_types_synonyms{SPLICESITE} = {'splicesite', 'splice_site', 'near_splice'}; 


genome_types{ENCODE} = 'encode'; % encode 1% regions (old encode) 
genome_types{GENOMIC} = 'genomic'; 


coding_allele_types = { ...
    {'all_synonymous', {'coding_synonymous', 'coding_synonymous_near_splice'}}, ...
    {'all_coding', {'coding_notMod3', 'coding_notMod3_near_splice', 'coding_synonymous', 'coding_synonymous_near_splice', ...
    'missense', 'missense_near_splice', 'splice_3', 'splice_5', 'stop_gained', 'stop_gained_near_splice', ...
    'stop_lost', 'stop_lost_near_splice'}}, ...
    {'all_missense', {'missense', 'missense_near_splice'}}, ... 
    {'all_nonsynonymous', {'coding_notMod3', 'coding_notMod3_near_splice', 'missense', 'missense_near_splice', ...
    'stop_gained', 'stop_gained_near_splice', 'stop_lost', 'stop_lost_near_splice'}}, ...
    {'all_lethal', {'coding_notMod3', 'coding_notMod3_near_splice',  'stop_gained', 'stop_gained_near_splice'}}, ...
    {'all_non_coding', {'intron', 'intergenic', 'near_gene_3', 'near_gene_5', 'utr_3', 'utr_5'}}, ...
    {'all_splicing',  {'coding_notMod3_near_splice', 'coding_synonymous_near_splice', 'missense_near_splice', ...
    'splice_3', 'splice_5', 'stop_gained_near_splice', 'stop_lost_near_splice'}} }; % Group different mutation classes together

USE_PI = 0; USE_OMEGA = 1; USE_ONLY_OMEGA = 2; % which conservation statistic to use

% stand
POS_STRAND = 1; REV_STRAND = 0; BOTH_STRANDS = 2; 
DNA_5TO3PRIME = 1; DNA_3TO5PRIME = 0;

% different organisms
ANCESTRAL=0;
HUMAN=1;

block_size = 1000000; % standard default block size
mid_block_size = 50; % block size for unions

% methods for ranking kmers
CONSERVATION = 0;
ENRICHMENT = 1;
CONSERVATION_AND_ENRICHMENT = 2;
CONSERVATION_TWO_CLADES = 3; 
CONSERVATION_AND_POLYMORPHISM = 4; 

USE_PI = 10; % here we use pi directly and not k-mers
KNOWN_PWMS = 20; % here simply search for enrichment of known transfac/other motifs
MOTIF_INSTANCE_PREDICTION = 100; % here we do not look for new motif but find instances of known motifs


% different regions types
ENRICHED_REGIONS = 0; % search for something in these regions vs. some background distribution
DISCRIMINATIVE_REGIONS = 1; % Here we want motifs that discriminate one set of regions from other sets of regions
RANKING_REGIONS = 2; % Here regions are given as a ranked list
SCORES_REGIONS = 3; % Here each region has some score associated with it (e.g. chip-seq data, expression etc.)

% Vector translating string of methods to/from numeric
motif_method_vec  = {CONSERVATION, 'CONSERVATION'; ...
    ENRICHMENT, 'ENRICHMENT'; ...
    CONSERVATION_AND_ENRICHMENT, 'CONSERVATION_AND_ENRICHMENT'; ...
    CONSERVATION_TWO_CLADES, 'CONSERVATION_TWO_CLADES'; ...
    CONSERVATION_AND_POLYMORPHISM, 'CONSERVATION_AND_POLYMORPHISM'; ...
    USE_PI, 'USE_PI'; ...
    KNOWN_PWMS, 'KNOWN_PWMS'; ...
    MOTIF_INSTANCE_PREDICTION, 'MOTIF_INSTANCE_PREDICTION'; ...
    ENRICHED_REGIONS, 'ENRICHED_REGIONS'; ... % here starting the regions types
    DISCRIMINATIVE_REGIONS, 'DISCRIMINATIVE_REGIONS'; ...
    RANKING_REGIONS, 'RANKING_REGIONS'; ...
    SCORES_REGIONS, 'SCORES_REGIONS'};

% Different ways of annotating genes
GENE_SYMBOL = 0;
REFSEQ = 1;
AFFYMETRIX = 2;

% Ordering of species in bit-strings given in different data files
POUYA_SPECIES_ORDER = {'ZEBRAFISH', 'MEDAKA', 'STICKLEBACK', 'FUGO', 'TETRAODON', 'FROG', 'LIZARD', ...
    'CHICKEN', 'PLATYPUS', 'MONODELPHIS', 'TENREC', 'ELEPHANT', 'ARMADILLO', 'COW', 'HORSE', 'CAT', 'DOG', ...
    'HEDGEHOG', 'SHREW', 'RABBIT', 'GUINEAPIG', 'MOUSE', 'RAT', 'TREESHREW', 'BUSHBABY', 'RHESUS', 'CHIMP'};

MANUEL_SPECIES_ORDER = {'HUMAN', 'CAVIA', 'TENREC', 'RABBIT', 'ELEPHANT', 'HEDGEHOG', 'SHREW', 'ARMADILLO', 'CAT', 'BAT', ...
    'SQUIRREL', 'PIKA', 'CHIMP', 'TREESHREW', 'BUSHBABY', 'DOG', 'RHESUS', 'MOUSE', 'RAT', 'MOUSELEMUR'};

% Note: Need to fill this with the correct one!!!!
NEW44_SPECIES_ORDER = { 'HUMAN', 'CHIMP', 'GORILLA', 'ORANGUTAN', 'RHESUS', 'MARMOSET', 'TARSIER', 'MOUSE_LEMUR', ...
    'BUSHBABY', 'TREESHREW', 'MOUSE', 'RAT', 'KANGAROO_RAT', 'GUINEA_PIG', 'SQUIRREL', 'RABBIT', ...
    'PIKA', 'ALPACA', 'DOLPHIN', 'COW', 'HORSE', 'CAT', 'DOG', 'MICROBAT', ...
    'MEGABAT', 'HEDGEHOG', 'SHREW', 'ELEPHANT', 'ROCK_HYRAX', 'TENREC', 'ARMADILLO', 'SLOTH', ...
    'MONODELPHIS', 'PLATYPUS', 'CHICKEN', 'ZEBRA_FINCH', 'LIZARD', 'X_TROPICALIS', 'TETRAODON', 'FUGU', ...
    'STICKLEBACK', 'MEDAKA', 'ZEBRAFISH', 'LAMPREY'};

% manuel_java_dir = '/seq/orzuk/24mammals/src/manuel_java/subversion/trunk/';  % directory with all the java code
manuel_java_dir = '/seq/mgarber/tools/'; % new dir with all java code

genesets_dir = '../data/gene_sets/'; % directory with all gene sets

% A vector of genome versions for each species. This agrees with manuels ordering (the ordering shouldn't matter here !!! )
% genome_versions_vec = { {'HUMAN', 'hg16', 'hg17', 'hg18'}, ...
%     {'CAVIA'}, ...
%     {'TENREC'}, ...
%     {'RABBIT'}, ...
%     {'ELEPHANT'}, ...
%     {'HEDGEHOG'}, ...
%     {'SHREW'}, ...
%     {'ARMADILLO'}, ...
%     {'CAT', 'felCat2'}, ...
%     {'BAT'}, ...
%     {'SQUIRRREL'}, ...
%     {'PIKA'}, ...
%     {'CHIMP', 'panTro1', 'panTro2'}, ...
%     {'TREESHREW'}, ...
%     {'BUSHBABY'}, ...
%     {'DOG', 'canFam1', 'canFam2'}, ...
%     {'RHESUS', 'rheMac2'}, ...
%     {'MOUSE', 'mm7', 'mm8', 'mm9'}, ...
%     {'RAT', 'rn2', 'rn3', 'rn4'}, ...
%     {'MOUSELEMUR'}};
%

% Format is: species genome-versions num-chroms which verterbrate, which clade (for human we don't want yet to use hg19 so it's not last)
% Megabat is not in encode! 
genome_versions_vec = {{'HUMAN', {'hg16', 'hg17', 'hg19', 'hg18', 'hg19'}, 23, {'MAMMAL', 'EUTHERIAN', 'PRIMATE', 'EUARCHONTOGLIRES', 'FULL', 'ENCODE'}}, ...
    {'CHIMP', {'panTro1', 'panTro2', 'panTro3', 'panTro4'}, 24, {'MAMMAL', 'EUTHERIAN', 'PRIMATE', 'EUARCHONTOGLIRES', 'FULL', 'ENCODE'}}, ...
    {'GORILLA', {'gorGor1', 'gorGor2', 'gorGor3'}, 24, {'MAMMAL', 'EUTHERIAN', 'PRIMATE', 'EUARCHONTOGLIRES', 'FULL', 'EXCLUDED', 'ENCODE'}}, ...
    {'BABOON', {'papHam1'}, 24, {'MAMMAL', 'EUTHERIAN', 'PRIMATE', 'EUARCHONTOGLIRES', 'FULL', 'EXCLUDED'}}, ... % not clear num. chroms !!! 
    {'ORANGUTAN', {'ponAbe2'}, 24, {'MAMMAL', 'EUTHERIAN', 'PRIMATE', 'EUARCHONTOGLIRES', '2X', 'EXCLUDED', 'ENCODE'}}, ...
    {'RHESUS', {'rheMac2', 'rheMac3'}, 24, {'MAMMAL', 'EUTHERIAN', 'PRIMATE', 'EUARCHONTOGLIRES', 'FULL', 'ENCODE'}}, ...
    {'MARMOSET', {'calJac1', 'calJac2', 'calJac3'}, -1, {'MAMMAL', 'EUTHERIAN', 'PRIMATE', 'EUARCHONTOGLIRES', '2X', 'EXCLUDED', 'ENCODE'}}, ...
    {'TARSIER', {'tarSyr1'}, -1, {'MAMMAL', 'EUTHERIAN', 'PRIMATE', 'EUARCHONTOGLIRES', '2X'}}, ...
    {'MOUSE_LEMUR', {'micMur1'}, -1, {'MAMMAL', 'EUTHERIAN', 'PRIMATE', 'EUARCHONTOGLIRES', '2X', 'ENCODE'}}, ...
    {'BUSHBABY', {'otoGar1', 'otoGar2', 'otoGar3'}, -1, {'MAMMAL', 'EUTHERIAN', 'PRIMATE', 'EUARCHONTOGLIRES', '2X'}}, ...
    {'TREESHREW', {'tupBel1'}, -1, {'MAMMAL', 'EUTHERIAN', 'EUARCHONTOGLIRES', '2X', 'ENCODE'}}, ...
    {'MOUSE', {'mm7', 'mm8', 'mm9', 'mm10'}, 20, {'MAMMAL', 'EUTHERIAN', 'RODENT', 'EUARCHONTOGLIRES', 'FULL', 'ENCODE'}}, ...
    {'RAT', {'rn2', 'rn3', 'rn4', 'rn5'}, 21, {'MAMMAL', 'EUTHERIAN', 'RODENT', 'EUARCHONTOGLIRES', 'FULL', 'ENCODE'}}, ...
    {'KANGAROO_RAT', {'dipOrd1'}, -1, {'MAMMAL', 'EUTHERIAN', 'RODENT', 'EUARCHONTOGLIRES', '2X'}}, ...
    {'GUINEA_PIG', {'cavPor3'}, -1, {'MAMMAL', 'EUTHERIAN', 'RODENT', 'EUARCHONTOGLIRES', 'FULL', 'ENCODE'}}, ...
    {'SQUIRREL', {'speTri1', 'speTri2'}, -1, {'MAMMAL', 'EUTHERIAN', 'RODENT', 'EUARCHONTOGLIRES', '2X', 'ENCODE'}}, ...
    {'RABBIT', {'oryCun1', 'oryCun2'}, 22, {'MAMMAL', 'EUTHERIAN', 'LAGOMORPH', 'EUARCHONTOGLIRES', '2X', 'ENCODE'}}, ...
    {'PIKA', {'ochPri2'}, -1, {'MAMMAL', 'EUTHERIAN', 'LAGOMORPH', 'EUARCHONTOGLIRES', '2X'}}, ...
    {'ALPACA', {'vicPac1'}, -1, {'MAMMAL', 'EUTHERIAN', 'LAURASIATHERIA', '2X'}}, ...
    {'DOLPHIN', {'turTru1', 'turTru2'}, -1, {'MAMMAL', 'EUTHERIAN', 'LAURASIATHERIA', '2X'}}, ...
    {'COW', {'bosTau4'}, 30, {'MAMMAL', 'EUTHERIAN', 'LAURASIATHERIA', 'FULL', 'ENCODE'}}, ...
    {'HORSE', {'equCab2'}, 32, {'MAMMAL', 'EUTHERIAN', 'LAURASIATHERIA', 'FULL', 'ENCODE'}}, ...
    {'CAT', {'felCat2', 'felCat3', 'felCat4', 'felCat5'}, 19, {'MAMMAL', 'EUTHERIAN', 'LAURASIATHERIA', 'CARNIVORA', '2X', 'ENCODE'}}, ...
    {'DOG', {'canFam1', 'canFam2', 'canFam3'}, 39, {'MAMMAL', 'EUTHERIAN', 'LAURASIATHERIA', 'CARNIVORA', 'FULL', 'ENCODE'}}, ...
    {'MICROBAT', {'myoLuc1', 'myoLuc2'}, -1, {'MAMMAL', 'EUTHERIAN', 'LAURASIATHERIA', 'CHRIOPTERA', '2X', 'ENCODE'}}, ...
    {'MEGABAT', {'pteVam1'}, -1, {'MAMMAL', 'EUTHERIAN', 'LAURASIATHERIA', 'CHRIOPTERA', '2X'}}, ...
    {'HEDGEHOG', {'eriEur1'}, 44, {'MAMMAL', 'EUTHERIAN', 'LAURASIATHERIA', '2X', 'ENCODE'}}, ...
    {'SHREW', {'sorAra1'}, -1, {'MAMMAL', 'EUTHERIAN', 'LAURASIATHERIA', '2X', 'ENCODE'}}, ...
    {'ELEPHANT', {'loxAfr2', 'loxAfr3'}, 28, {'MAMMAL', 'EUTHERIAN', 'AFROTHERIA', 'ATLANTOGENATA', '2X', 'ENCODE'}}, ...
    {'ROCK_HYRAX', {'proCap1'}, -1, {'MAMMAL', 'EUTHERIAN', 'AFROTHERIA', 'ATLANTOGENATA', '2X', 'ENCODE'}}, ...
    {'TENREC', {'echTel1'}, -1, {'MAMMAL', 'EUTHERIAN', 'AFROTHERIA', 'ATLANTOGENATA', '2X', 'ENCODE'}}, ...
    {'ARMADILLO', {'dasNov2'}, -1, {'MAMMAL', 'EUTHERIAN', 'EDENTATA', 'ATLANTOGENATA', '2X', 'ENCODE'}}, ...
    {'SLOTH', {'choHof1'}, -1, {'MAMMAL', 'EUTHERIAN', 'EDENTATA', 'ATLANTOGENATA', '2X', 'ENCODE'}}, ...
    {'PANDA', {'ailMel'}, -1, {'MAMMAL'}}, ... % New mammals from 60way alignment
    {'SHEEP', {'oviAri1'}, -1, {'MAMMAL'}}, ... % New mammals from 60way alignment
    {'SQUIRREL_MONKEY', {'saiBol1'}, -1, {'MAMMAL'}}, ... % New mammals from 60way alignment
    {'NAKED_MOLERAT', {'hetGla2'}, -1, {'MAMMAL'}}, ... % New mammals from 60way alignment
    {'GIBBON', {'nomLeu2'}, -1, {'MAMMAL', 'PRIMATE'}}, ... % New mammals from 60way alignment
    {'TURTLE', {'chrPic1'}, -1, {'REPTILE'}}, ... % New reptile from 60way alignment
    {'COD', {'gadMor1'}, -1, {'FISH'}}, ... % New fish from 60way alignment
    {'COELANCANTH', {'latCha1'}, -1, {'FISH'}}, ... % New fish from 60way alignment    
    {'TILAPIA',  {'oreNil2'}, -1, {'FISH'}}, ... % New fish from 60way alignment    
    {'MONODELPHIS', {'monDom4', 'monDom5'}, 11, {'MAMMAL', 'MARSUPIAL', 'OPOSSUM', 'FULL', 'ENCODE'}}, ...
    {'PLATYPUS', {'ornAna1'}, 26, {'MAMMAL', 'MONOTREME', 'FULL', 'ENCODE'}}, ...
    {'CHICKEN', {'galGal3'}, 39, {'BIRD', 'FULL', 'ENCODE'}}, ...
    {'ZEBRA_FINCH', {'taeGut1'}, -1, {'BIRD', 'FULL', 'ENCODE'}}, ...
    {'LIZARD', {'anoCar1', 'anoCar2'}, -1, {'REPTILE', 'FULL', 'ENCODE'}}, ...
    {'X_TROPICALIS', {'xenTro2', 'xenTro3'}, -1, {'AMPHIBIAN', '2X', 'ENCODE'}}, ...
    {'TETRAODON', {'tetNig1', 'tetNig2'}, -1, {'FISH', 'FULL', 'ENCODE'}}, ...
    {'FUGU', {'fr2', 'fr3'}, -1, {'FISH', 'FULL', 'ENCODE'}}, ...
    {'STICKLEBACK', {'gasAcu1'}, -1, {'FISH', 'FULL', 'ENCODE'}}, ...
    {'MEDAKA', {'oryLat2'}, -1, {'FISH', '2X', 'ENCODE'}}, ...
    {'ZEBRAFISH', {'danRer5', 'danRer6', 'danRer7'}, -1, {'FISH', 'FULL', 'ENCODE'}}, ...
    {'LAMPREY', {'petMar1'}, -1, {'FISH', '2X', 'EXCLUDED', 'ENCODE'}}};

genome_versions_encode_vec = cell(length(genome_versions_vec), 1); % New! add encode names for species! they are differnet sometimes from the 'standard' name! 
for i_org_str=1:length(genome_versions_vec)
    genome_versions_encode_vec{i_org_str} = genome_versions_vec{i_org_str}{1};    
end
genome_versions_encode_vec{strmatch('SQUIRREL', genome_versions_encode_vec, 'exact')} = 'ST_SQUIRREL'; % enumerate all differences between encode and 2x names
genome_versions_encode_vec{strmatch('RHESUS', genome_versions_encode_vec)}  = 'MACAQUE'; % enumerate all differences between encode and 2x names
genome_versions_encode_vec{strmatch('BUSHBABY', genome_versions_encode_vec)}  = 'BUSH_BABY'; % enumerate all differences between encode and 2x names
genome_versions_encode_vec{strmatch('TREESHREW', genome_versions_encode_vec)}  = 'TREE_SHREW'; % enumerate all differences between encode and 2x names
genome_versions_encode_vec{strmatch('MICROBAT', genome_versions_encode_vec)}  = 'SBBAT'; % enumerate all differences between encode and 2x names
% genome_versions_encode_vec{strmatch('MEGABAT', genome_versions_encode_vec)}  = 'RFBAT'; % enumerate all differences between encode and 2x names - this is a different bat! we exclude it!!! 



% TREE: (((((((((((((human:0.008662,chimp:0.008831):0.032232,macaque:0.046182):0.112209,tarsier:0.146881):0.013022,(mouse_lemur:0.1209
% 39,bush_baby:0.153638):0.039735):0.015825,tree_shrew:0.206994):0.005326,(((((mouse:0.084059,rat:0.086027):0.230338,kangaroo_rat:0.26
% 6761):0.020554,guinea_pig:0.236908):0.007059,st_squirrel:0.186292):0.036109,(rabbit:0.157475,pika:0.236918):0.123221):0.015574):0.02
% 6874,(((alpaca:0.122102,(dolphin:0.079064,cow:0.155441):0.030294):0.044851,((horse:0.133905,(cat:0.122675,dog:0.123578):0.062649):0.
% 003921,(sbbat:0.170326,rfbat:0.148669):0.051961):0.003673):0.014867,(hedgehog:0.275009,shrew:0.391209):0.04875):0.023653):0.026406,(
% ((elephant:0.107652,rock_hyrax:0.197456):0.038652,tenrec:0.322578):0.065335,(armadillo:0.132814,sloth:0.120099):0.070272):0.010425):
% 0.319842,monDom4:0.482272):0.066275,ornAna1:0.571394):0.130625,((galGal3:0.240715,taeGut1:0.21643):0.317015,anoCar1:0.984317):0.1674
% 39):0.168382,xenTro2:1.765299):0.37068,(((tetNig1:0.269914,fr2:0.24186):0.259502,(gasAcu1:0.391367,oryLat2:0.645255):0.099892):0.542
% 691,danRer5:1.364005):0.34791):0.1


excluded_species = {'GORILLA', 'MARMOSET', 'ORANGUTAN', 'LAMPREY'}; % These are not in the ucsc tree and 2x paper (also Lamprey)
% More genomes we want (we don't have them in the alignment, but perhaps in the future):
% sequenced: PANDA, KANGAROO, PIG (Porcine, partial), TILAPIA (partial)
% not: MOLERAT, SHEEP, BABOON, WHALE, CAMEL, TURTULE, SKATE (dag
% dalton), SPOTTED_GAR (fish), PEROMYSCUS (american deer mice), BONOBO
% (like-chimp), ELEPHANT_SHREW, FERRET, GIBBON, FLYING_LEMUR, LLAMA
% (cancelled), CRAB_EATING_MACAQUE, MOLE (on hold), PANGOLIN (on hold),
% PRAIRIE_VOLE, SQUIRREL_MONKEY, VERVET (monkey), WALLABY
% GIRAFFE, PORCUPINE, BEAVER, SEAL, BEAR, SEA_LION, MONGOOSE, FIREFOX,
% SKUNK, RACOON, WEASEL, WOLF, FOX, COYOTA, JACKAL, HYENA, LION, TIGER,
% JAGUAR, LEOPARD, VIVERRIDAE, PORPOISE, TAPIR, HIPPO, DEER, RHINO, GOAT,
% SOLENODON, ANTEATER (David Haussler plans to do almost ALL of them - a
% few thousands)

% Num chroms (including one for the sex chroms). This matches manuels species names. -1 means we don't know yet
NUM_SPECIES_CHROMS = [23, -1, -1, 22, -1, -1, -1, -1, 19, -1, ...
    -1, -1, 24, -1, -1, 39, 24, 20, 21, -1];

% known_pwms_file = '../data/pwms_small.mat'; % file containing all known pwms: transfac, xie's, jaspar(?) etc.
% known_pwms_file = '../data/pwms_standard.mat'; % file containig known pwms we can use
% known_pwms_file = '../data/pwms_clustered.mat'; % file containing clustered known pwms
% known_pwms_file = '../data/pwms_union.mat'; % file containing all known pwms: transfac, xie's, jaspar(?) etc.
% known_pwms_file = '../data/pwms/pwms_noa.mat';
if(~exist('known_pwms_file', 'var') || isempty(known_pwms_file))
    known_pwms_file = '../data/pwms_new_union_set.mat';
end
% known_pwms_file = '../data/pwms/pwms_irf.mat';
derich_alpha = 0.05; % default derichlet correction for pwms etc.

% Parameters of the evolutionary model,
PI_0 = [0.3 0.2 0.2 0.3]; % AT rich
Q_0 = [ -1.0860163084198684 0.25225991842162027  0.6937865599085344 0.13996983008971356; ...
    0.16817327894774686 -0.9225558555676914  0.2866374459445824 0.46774513067536205;...
    0.4625243732723563  0.2866374459445824 -0.9192377889311856 0.17007596971424688;...
    0.13996983008971356   0.701617696013043  0.2551139545713703 -1.0967014806741269];% rate matrix


PHAST=1; PHYLOP=2; SIPHY=3; GERP=4;  % Methods for conservation


