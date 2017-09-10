% Here get the height of each column to be according to its information
% content
function string1 = draw_pwm_new(pwm, add_legend, varargin)

L = length(pwm);
string1 = cell(1,L);
if(~exist('add_legend', 'var'))
    add_legend = 1;
end


% temp = max(pwm, [], 1);
% A = 1; 
% C = 2; 
% G = 3; 
% T = 4; 
% indexes = zeros(1,L);
% for i = 1:L
%     indy = find(pwm(:,i) == temp(i));
%     indy = indy(1);
%     indexes(i) = indy;
% 
%     switch indy
%         case 1
%             string1{i} = 'A';
%         case 2
%             string1{i} = 'C';
%         case 3
%             string1{i} = 'G';
%         case 4
%             string1{i} = 'T';
%     end
% end

string1 = mat2cell(pwm_to_consensus(pwm), 1, ones(1,L));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%% New way  %%%%%%%%%%%%%%%%%%%%%%%%%%
heights = pwm;
heights(3,:) = heights(3,:) + heights(4,:);
heights(2,:) = heights(2,:) + heights(3,:);
heights(1,:) = 2-entropy(pwm);
for i=2:4
    heights(i,:) = heights(i,:) .* heights(1,:);
end

%figure;
hold on; bar(heights(1,:)); bar(heights(2,:), 'y'); bar(heights(3,:), 'g'); bar(heights(4,:), 'r');
if(add_legend)
    legend('A', 'C', 'G', 'T'); % xlabel('consensus');
end
ylabel('Information'); title('Position weight Matrix');
set(gca,'xtick',[1:L],'xticklabel', string1);



dra = 1;