% Clipp gaussian variables into discrete variables 
function CLIPPED_VEC = ClipModelObservedData( HMM_MODEL, Y_VEC)

Threshold = TwoGaussiansThreshold( HMM_MODEL.MEW(1), HMM_MODEL.SIGMA(1),HMM_MODEL.MEW(2),HMM_MODEL.SIGMA(2));
CLIPPED_VEC = (Y_VEC > Threshold);

