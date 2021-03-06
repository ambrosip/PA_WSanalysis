function [lightEvokedCurrents,dataPerSweepCh1,dataPerSweepCh2] = mono(obj, varargin)

[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);

% optional arguments
numvarargs = length(varargin);
optargs = {2 0.005 1.94 2.14 -3000 500};
optargs(1:numvarargs) = varargin;
[lightOnsetTime, lightDuration, xmin, xmax, ymin, ymax] = optargs{:};

totalActiveChannels = 2;
sampleRate = obj.header.Acquisition.SampleRate;
lightOnsetDataPoint = lightOnsetTime*sampleRate;
lightOffDataPoint = (lightOnsetTime+lightDuration)*sampleRate;
afterLightDataPoint = (lightOnsetTime+0.1)*sampleRate;

lightEvokedCurrents = [];
dataPerSweepCh1 = [];
dataPerSweepCh2 = [];

% checking for incomplete sweeps and not analyzing incomplete sweeps (to
% avoid this error: "Dimensions of matrices being concatenated are not
% consistent." when calculating this: "dataPerSweepCh2 = [dataPerSweepCh2,
% ych2]" and other concatenations.
if numel(fieldnames(obj.sweeps)) < obj.header.NSweepsPerRun
    lastSweepNumber = lastSweepNumber - (obj.header.NSweepsPerRun - numel(fieldnames(obj.sweeps)));
    allSweeps = firstSweepNumber:lastSweepNumber;
end 
    
    % plotting one figure with all sweeps and average data 
    figure('name', strcat(obj.file,' (all) - 2 channels')); % naming figure file
    
    % code for channel 2
    for channel = 2:totalActiveChannels
        for sweepNumber = allSweeps
        hold on;
        [xch2,ych2] = obj.xy(sweepNumber, channel);     
        subplot(totalActiveChannels,1,channel)
        plot(xch2,ych2,'k');
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([xmin xmax yminhere ymaxhere])
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(channel), ' (', obj.header.Acquisition.AnalogChannelUnits(channel), ')'));
        dataPerSweepCh2 = [dataPerSweepCh2, ych2];
        end
        xlabel('Time (s)');
        hold off;
    end
    
    % separate code for channel 1 to apply axis range and figure title
    for sweepNumber = allSweeps  
        hold on;
        [x,y] = obj.xy(sweepNumber, 1);
        subplot(totalActiveChannels,1,1)
        plot(x,y,'Color', [0.75, 0.75, 0.75]);
        axis([xmin xmax ymin ymax])
%         xlabel('Time (s)');
        ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
        title([obj.file ' (all) - 2 channels'],'Interpreter','none');
        
        dataPerSweepCh1 = [dataPerSweepCh1, y];        
        lightEvokedCurrents = [lightEvokedCurrents, min(y(lightOnsetDataPoint:afterLightDataPoint))];        
    end

    % plotting average for channel 1
    % dataPerSweepCh1 has a column for each sweep. mean(dataPerSweepCh1,2)
    % calculates the average column. mean(dataPerSweepCh1) calculates the
    % wrong thing: the average row.
    subplot(totalActiveChannels,1,1)
    plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
    rectangle('Position', [lightOnsetTime 0 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
%     line(lightOnsetDataPoint:lightOffDataPoint,100:100,'Color','blue','LineWidth',3)
    axis([xmin xmax ymin ymax])
    hold off;    
    
    movegui('northwest');
    
    % plotting figure with only Ch1
    figure('name', strcat(obj.file,' (all)')); % naming figure file
    hold on;
    plot(x, dataPerSweepCh1,'Color', [0.75, 0.75, 0.75]);
    plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
    rectangle('Position', [lightOnsetTime 0 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
    axis([xmin xmax ymin ymax])
    xlabel('Time (s)');
    ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
    title([obj.file ' (all)'],'Interpreter','none');
    hold off;
    
    
    % plotting figure with only Ch1 niceplot
    figure('name', strcat(obj.file,' (all) - niceplot')); % naming figure file
    hold on;
    plot(x, dataPerSweepCh1,'Color', [0.75, 0.75, 0.75]);
    plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
    rectangle('Position', [lightOnsetTime 200 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
    axis([lightOnsetTime-0.01 lightOnsetTime+0.1 -2000 500]) % used to be -7000
    xlabel('Time (s)');
    ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
    title([obj.file ' (all)'],'Interpreter','none');
    hold off;
%     set(gca,'Visible','off')
    

end