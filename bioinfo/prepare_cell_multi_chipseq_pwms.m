% Prepare some pwms which are not present in our collection but appear to
% show up in the cell paper of 13 TFs of the Singapore group 
KLF4 = iupac_to_pwm('GGGTGTGCC'); 
KLF4(:,4) = [0 0.3 0 0.7]';
KLF4(:,6) = [0 0 0.33 0.67]';
KLF4(:,8) = [0 0 0.9 0.1]';
KLF4(:,9) = [0.05 0.7 0 0.25]';
KLF4(:,10) = [0 0.6 0.05 0.35]';
seqlogo(KLF4);

ESRRB = iupac_to_pwm('GGTCAAGGTCAC');
ESRRB(:,1) = [0.13 0.05 0.7 0.16]';
ESRRB(:,2) = [0.1 0.2 0.54 0.16]';
ESRRB(:,3) = [0 0.2 0 0.8]';
ESRRB(:,4) = [0 0.85 0.15 0]';
ESRRB(:,end) = [0.2 0.4 0.2 0.2]';
seqlogo(ESRRB);

ZFX = iupac_to_pwm('CGCNAGGCCGCG');
ZFX(:,1) = [0.15 0.55 0.25 0.05]';
ZFX(:,2) = [0.05 0.18 0.7 0.07]';
ZFX(:,3) = [0 0.9 0.1 0]';
ZFX(:,4) = [0.2 0.29 0.21 0.3]';
ZFX(:,5) = [0.55 0.15 0.2 0.1]';
ZFX(:,10) = [0.04 0.15 0.65 0.16]';
ZFX(:,11) = [0.07 0.7 0.18 0.05]';
ZFX(:,12) = [0.07 0.05 0.8 0.08]';

seqlogo(ZFX);

load('../data/pwms_union.mat');
n = size(pwms, 1)
pwms{end+1,1} = 'KLF4'; pwms{end,2} = KLF4; pwms{end,3} = ''; pwms{end,4} = '';
pwms{end+1,1} = 'ESRRB'; pwms{end,2} = ESRRB; pwms{end,3} = ''; pwms{end,4} = '';
pwms{end+1,1} = 'ZFX'; pwms{end,2} = ZFX; pwms{end,3} = ''; pwms{end,4} = '';
% save('../data/pwms_union.mat', 'pwms'); % Do it only once !!! 


