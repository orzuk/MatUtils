global BNT_HOME
%BNT_HOME = 'C:\kpmurphy\matlab\FullBNT';
% BNT_HOME = '/home/ai2/murphyk/matlab/FullBNT';
BNT_HOME = 'E:\networks\BN software\Kevin Murphy\FullBNT';
v = version;
if v(1)=='5'
  addpath(genpath(BNT_HOME,0))
else
  addpath(genpath(BNT_HOME))
end
% fails to add directories which only contain directories but no regular files
% e.g., BNT/inference. Hence we add dummy files to such directories.
% This bug has been fixed in matlab 6.5

