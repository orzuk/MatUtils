% Generate a pwm from genomic sites
% 
% The input: 
% str - strings representing genomic sequences
% derich - derichlet correction 
% letters_flag - whether input is alphabetic or numeric
% 
% The output: 
% pwm - the pwm 
% bs - the binding sites
% L - pwm length
% 
function [pwm, bs, L] = sites_to_pwm(str, derich, letters_flag)


% fin=str;
% nLines=0;
% string = [];
% 
% i=0;
% while (1==1) %nLines < 50000)
%     line=fgetl(fin);
%     if (~ischar(line)) %EOF
%         break;
%     end
%     
%     % Set L according to the first line
%     if(i==0)
%         L = length(line);
%     end 
%     i=i+1;
%     
%     nLines=nLines+1;
%     if mod(nLines,2000)==0
%         disp(nLines);
%     end
%     if (~(isempty(line)) & (~(strcmp(line(1),'<'))))
%         string = [string line];
%     end
% %     tabs=find(line==9); %9 is TAB
% %     starts=[1 tabs+1];
% %     ends=[tabs-1 length(line)];
% % 
% %     for i=1:length(starts)
% %         table{nLines, i}=line(starts(i):ends(i));
% %     end
% end
% 
% 
% fclose(fin);



% Transfer sequence to numbers ..


if(letters_flag)
    genome = double(str);


    genome(find(genome == double('a'))) = 1;
    genome(find(genome == double('t'))) = 2;
    genome(find(genome == double('c'))) = 3;
    genome(find(genome == double('g'))) = 4;
    genome(find(genome == double('A'))) = 1;
    genome(find(genome == double('T'))) = 2;
    genome(find(genome == double('C'))) = 3;
    genome(find(genome == double('G'))) = 4;
else
    genome = str;
end
    

% Create a matrix
% num_sites = length(genome)/L;
% bs = reshape(genome,  L, num_sites);
% bs = bs'; 
num_sites = size(genome,1);
L=size(genome,2);
bs = genome;

pwm = zeros(4, L);

for i=1:L
    for j=1:4
        pwm(j, i) = length( find( bs(:,i) == j ) );
    end
end

pwm;

% A relative derichle correction : 
pwm = pwm / num_sites;
pwm = (pwm + derich )/ ( 1 + 4*derich);


 



