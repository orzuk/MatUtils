% Plot various difference metrics between two HMMs 
function [DUMMY] = PlotModelsDiff( MODELS, TRUE_MODEL)


SS_M = zeros(1, length(MODELS));
SS_MEW = zeros(1, length(MODELS));
SS_SIGMA = zeros(1, length(MODELS));
for i=1:length(MODELS)
    
    
    % Compute the sum of squares of differences between models 
    SS_M(i) = sum(sum((MODELS{i}.M - TRUE_MODEL.M).^2));
    SS_MEW(i) = sum(sum((MODELS{i}.MEW - TRUE_MODEL.MEW).^2));
    SS_SIGMA(i) = sum(sum((MODELS{i}.SIGMA - TRUE_MODEL.SIGMA).^2));
end


SS_M = SS_M / length(TRUE_MODEL.PI)*(length(TRUE_MODEL.PI)-1);
SS_SIGMA = SS_SIGMA / length(TRUE_MODEL.PI);
SS_MEW = SS_MEW / length(TRUE_MODEL.PI); 

figure; plot([1:length(MODELS)], SS_M, '.'); title('M Sum Of Squares Error');  xlabel('Model'); ylabel('Error'); 
figure; plot([1:length(MODELS)], SS_MEW, '.'); title('MEW Sum Of Squares Error');  xlabel('Model'); ylabel('Error'); 
figure; plot([1:length(MODELS)], SS_SIGMA, '.'); title('SIGMA Sum Of Squares Error');  xlabel('Model'); ylabel('Error'); 

DUMMY = 0; 