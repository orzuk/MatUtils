% function AssignGeneralConstants()
% General constants to use 
beep off; % avoid annoying beeping 
UNIX=0; PC=1;% Check if machine is unix or pc
LEFT = 0; RIGHT = 1; % directions 
COLUMN = 0; ROW = 1;  % vector orientations 
BFS = 0; DFS = 1; % breadth-first-search vs. depth-first-search
epsilon = 0.00000000001;
BIG_NUM = 9999999999.9;
color_vec = 'brgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmcbrgkmc';
matlab_root = '/broad/tools/apps/matlab77/bin/matlab'; % matlab root dir 
RED = 1; GREEN = 2; BLUE = 3; MAGENTA = 4; ORANGE = 5; BLACK = 6; CYAN = 7; BROWN = 8; YELLOW = 9;
tab = sprintf('\t'); eoln = sprintf('\n'); 
MAT = 1; TXT = 2; MAT_AND_TXT = 3; % file formats 
BLOCK_SIZE = 1000000; % beyond that split intervals/arrays 
EMPTY = 1; CORRUPT = 2;  GOOD = 0; % file states 

AND = 1000; OR = 1001; XOR = 1002; NOT = 1003; SUM = 1004; EXP = 1005; % Gates
SIGMOID = 2000; ADDITIVE = 2001; LOGISTIC = 2002; MULTIPLICATIVE = 2003; 
AFFINE = 2004; 
MAJORITY = 3000; K_OF_N = 3001; K_OR_MORE_OF_N = 3002; 
precision = 3; % precision when displaying real numbers 


machine = ispc;
if(machine == UNIX)
    format_fig_vec = {'fig'};
else
    format_fig_vec = {'fig', 'epsc', 'jpg'}; % how to save matlab figures
end

