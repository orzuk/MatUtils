% Plot a mixture model and see how does it look !!!! 
function [DUMMY] = PlotModel(full_path_data_file_name, full_path_model_and_output_file_name, ChromosomeToPlot, KMin, num_points)

% First simulate the sequence 
OUT_VEC = SimulateSequenceFromModelMatlab(full_path_data_file_name, full_path_model_and_output_file_name, ChromosomeToPlot, KMin, num_points); 


% Now Plot the random points 
figure; hist(OUT_VEC, 100); title(' Model Simulations Histogram'); xlabel(' Expression '); ylabel(' Frequency'); 

figure; plot([1:num_points], OUT_VEC,  'm*'); title([' Model Simulations : Mean ' num2str(mean(OUT_VEC)) ' std. dv. ' num2str(std(OUT_VEC))]); 
    xlabel(' Location'); ylabel(' Frequency'); 

XMIN = 0; 
XMAX = num_points;
YMIN = -2; 
YMAX = KMin;

AXIS([XMIN XMAX YMIN YMAX]);

H = line([0, XMAX], [0, 0]); % set the boundary between p and q armsset
set(H, 'LineWidth', 2);  set(H, 'Color', 'k');   % Color it in black


DUMMY = 0; 