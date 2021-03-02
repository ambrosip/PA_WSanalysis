function [mean_data, preLightHzMean, duringLightHzMean, postLightHzMean] = firingSpaghettiPrep(mean_sd_data)

preLightHzMean = mean_sd_data(:,1);
duringLightHzMean = mean_sd_data(:,3);
postLightHzMean = mean_sd_data(:,5);
mean_data = [preLightHzMean duringLightHzMean postLightHzMean];

end