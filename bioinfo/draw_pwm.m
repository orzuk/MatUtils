% Plot a pwm. Similar to Matlab's seqlogo but gives a 'normal' plot window
% which can be incorporated with other functions, e.g. subplots)
function string1 = draw_pwm(pwm)

L = length(pwm);

temp = max(pwm, [], 1);

string1 = cell(1,L);

A = 1; 
C = 2; 
G = 3; 
T = 4; 

indexes = zeros(1,L);
for i = 1:L
    indy = find(pwm(:,i) == temp(i));
    indy = indy(1);
    indexes(i) = indy;

    switch indy
        case 1
            string1{i} = 'A';
        case 2
            string1{i} = 'C';
        case 3
            string1{i} = 'G';
        case 4
            string1{i} = 'T';
    end

end



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
% subplot(3,3,i_ind);
hold on; bar('v6',heights(1,:)); bar('v6',heights(2,:), 'y'); bar('v6',heights(3,:), 'g'); bar('v6',heights(4,:), 'r');
legend('A', 'C', 'G', 'T'); % xlabel('consensus');
ylabel('Information'); title('Position weight Matrix');
set(gca,'xtick',[1:L],'xticklabel', string1);



dra = 1;