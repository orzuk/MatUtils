% function AssignGeneralConstants()
% General constants to use
beep off; % avoid annoying beeping
UNIX=0; PC=1;% Check if machine is unix or pc
LEFT = 0; RIGHT = 1; % directions
COLUMN = 0; ROW = 1;  % vector orientations
OLD = 0; NEW = 1;
BFS = 0; DFS = 1; % breadth-first-search vs. depth-first-search
epsilon = 0.00000000001;
BIG_NUM = 9999999999.9;
color_vec = 'brgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmc';
% symbol_vec = '.:-xo*s+dv^ph.:-xo*s+dv^ph.:-xo*s+dv^ph.:-xo*s+dv^ph.:-xo*s+dv^ph.:-xo*s+dv^ph.:-xo*s+dv^ph';
symbol_vec = {'-', ':', '--', '-.', 'o', 'x', 's', 'v', 'd', '*', '^', '.'}; 
symbol_str_vec = {'solid', 'dotted', 'dashed', '-.',  'o', 'square', 'v', 'diamond', 'star', '^', 'dot'}; 

%matlab_root = '/broad/tools/apps/matlab77/bin/matlab'; % matlab root dir
matlab_root = '/broad/tools/apps/matlab2009b/bin/matlab'; % matlab root dir (new version)
%matlab_root = '/broad/tools/apps/matlab2010a/bin/matlab'; % matlab root dir (new version)
RED = 1; GREEN = 2; BLUE = 3; MAGENTA = 4; ORANGE = 5; BLACK = 6; CYAN = 7; BROWN = 8; YELLOW = 9;
tab = sprintf('\t'); eoln = sprintf('\n');
MAT = 1; TXT = 2; MAT_AND_TXT = 3; % file formats
BLOCK_SIZE = 1000000; % beyond that split intervals/arrays
EMPTY = 1; CORRUPT = 2;  GOOD = 0; % file states

AND = 1000; OR = 1001; XOR = 1002; NOT = 1003; SUM = 1004;
MAX = 1005; MIN = 1006; EXP = 1007; THRESHOLD = 1008; RANDOM = 1009; % Gates
SIGMOID = 2000; ADDITIVE = 2001; LOGISTIC = 2002; MULTIPLICATIVE = 2003; AFFINE = 2004; LIABILITY = 2005; PROBIT = 2006;
LOGIT = LOGISTIC; % synonimous names
MAJORITY = 3000; K_OF_N = 3001; K_OR_MORE_OF_N = 3002;
DOMINANT = 4000; RECESSIVE = 4001;
precision = 3; % precision when displaying real numbers

MALE = 1; FEMALE = 0; % genders

machine = ispc;
%if(machine == UNIX)
%    format_fig_vec = {'fig'};
%else
    format_fig_vec = {'fig', 'epsc', 'jpg', 'pdf'}; % how to save matlab figures
%end

switch machine
    case UNIX
        matlab_libs_root_dir = '~orzuk/public_html/matlab/libs';
        matlab_word_size = 64; % NEW! assume we're working in 64 bit
    otherwise
        matlab_libs_root_dir = 'C:\research\matlab\libs'; % new laptop
        matlab_word_size = 32; % assume we're working in 32 bit
end
%                        matlab_libs_root_dir = 'Y:\public_html\matlab\libs';
