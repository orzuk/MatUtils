% Add a right y-axis with possibly different scale
%
% Input:
% right_y_scale - ticks for right y axis
% right_y_label_str - a label for right y axis
%
function add_right_yticks(right_y_scale, right_y_label_str)

%plot(rand(1,10));       %# Plot some random data
%ylabel(gca,'scale 1');  %# Add a label to the left y axis
set(gca,'Box','off');   %# Turn off the box surrounding the whole axes
axesPosition = get(gca,'Position');          %# Get the current axes position
hNewAxes = axes('Position',axesPosition,...  %# Place a new axes on top...
    'Color','none',...           %#   ... with no background color
    'YLim',right_y_scale,...     %#   ... and a different scale
    'YAxisLocation','right',...  %#   ... located on the right
    'XTick',[],...               %#   ... with no x tick marks
    'Box','off');                %#   ... and no surrounding box
if(exist('right_y_label_str', 'var') && (~isempty(right_y_label_str)))
    ylabel(hNewAxes, right_y_label_str);  %# Add a label to the right y axis
end
