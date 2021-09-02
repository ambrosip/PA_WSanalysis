function [preLightHzMean, preLightHzSD, duringLightHzMean, duringLightHzSD] = firingScatterPrep(data)

preLightHzMean = data(:,1);
preLightHzSD = 2*data(:,2);
duringLightHzMean = data(:,3);
duringLightHzSD = 2*data(:,4);

end
