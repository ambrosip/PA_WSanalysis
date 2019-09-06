function [lightEvokedCurrents,dataPerSweepCh1,dataPerSweepCh2,seriesResistance] = normmono(obj, varargin)

[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);

% optional arguments
numvarargs = length(varargin);
% optargs = {1.99 2 0.005 1.94 2.14 -3000 500 1 'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries\From MATLAB'};
optargs = {1.99 2 0.005 1.94 2.14 -3000 500 1 'D:\Temp\From MATLAB'};
% optargs = {1.99 2 0.005 1.94 2.14 -3000 500 1 'D:\Temp\Data summaries\2019-07-20 th flp naive'};
optargs(1:numvarargs) = varargin;
[baselineOnsetTime, lightOnsetTime, lightDuration, xmin, xmax, ymin, ymax, rsTestPulseOnsetTime, savefileto] = optargs{:};

% variables that will be used
totalActiveChannels = 2;
sampleRate = obj.header.Acquisition.SampleRate;
lightOnsetDataPoint = lightOnsetTime*sampleRate;
lightOffDataPoint = (lightOnsetTime+lightDuration)*sampleRate;
afterLightDataPoint = (lightOnsetTime+0.1)*sampleRate;
baselineOnsetDataPoint = baselineOnsetTime*sampleRate;
rsBaselineDataPointInterval = ((rsTestPulseOnsetTime-0.1)*sampleRate):(rsTestPulseOnsetTime*sampleRate);
rsFirstTransientDataPointInterval = (rsTestPulseOnsetTime*sampleRate):(rsTestPulseOnsetTime+0.0025)*sampleRate;

% matrixes that will be filled
lightEvokedCurrents = [];
dataPerSweepCh1 = [];
dataPerSweepCh2 = [];
seriesResistance = [];
data = [];
allRs = [];

    % checking for incomplete sweeps and not analyzing incomplete sweeps (to
    % avoid this error: "Dimensions of matrices being concatenated are not
    % consistent." when calculating this: "dataPerSweepCh2 = [dataPerSweepCh2,
    % ych2]" and other concatenations.
    if numel(fieldnames(obj.sweeps)) < obj.header.NSweepsPerRun
        lastSweepNumber = lastSweepNumber - (obj.header.NSweepsPerRun - numel(fieldnames(obj.sweeps)));
        allSweeps = firstSweepNumber:lastSweepNumber;
    end    
    
    % plotting one figure with all sweeps and average data 
    figure('name', strcat(obj.file,' (all) - 2 channels - Normalized')); % naming figure file
    
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
        
        % calculating series resistance
        rsBaselineCurrent = mean(y(rsBaselineDataPointInterval));
        rsTransientCurrent = min(y(rsFirstTransientDataPointInterval));
        dCurrent = rsTransientCurrent-rsBaselineCurrent;
        dVoltage = -5;
        seriesResistance = 1000*dVoltage/dCurrent; %mV/pA equals Gohm
        
        % normalizing data based on mean baseline current between t=1.99 (or whatever was the input time) and t=2
        y = y-mean(y(baselineOnsetDataPoint:lightOnsetDataPoint));
        
        subplot(totalActiveChannels,1,1)
        plot(x,y,'Color', [0.75, 0.75, 0.75]);
        axis([xmin xmax ymin ymax])
%         xlabel('Time (s)');
        ylabel(strcat("Normalized ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
        title([obj.file ' (all) - 2 channels - Normalized'],'Interpreter','none');
        
        dataPerSweepCh1 = [dataPerSweepCh1, y];
        lightEvokedCurrents = [lightEvokedCurrents, min(y(lightOnsetDataPoint:afterLightDataPoint))];  
        lightEvokedCurrent = min(y(lightOnsetDataPoint:afterLightDataPoint));
        
        data = [data; sweepNumber, seriesResistance, lightEvokedCurrent];
        allRs = [allRs, seriesResistance];
        
    end

    % plotting average for channel 1
    % dataPerSweepCh1 has a column for each sweep. mean(dataPerSweepCh1,2)
    % calculates the average column. mean(dataPerSweepCh1) calculates the
    % wrong thing: the average row.
    subplot(totalActiveChannels,1,1)
    plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
    rectangle('Position', [lightOnsetTime 200 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
%     line(lightOnsetDataPoint:lightOffDataPoint,100:100,'Color','blue','LineWidth',3)
    axis([xmin xmax ymin ymax])
    hold off;    
    
    movegui('northwest');
    
%     % plotting figure with only Ch1
%         figure('name', strcat(obj.file,' (all) - Normalized')); % naming figure file
%         hold on;
%         plot(x, dataPerSweepCh1,'Color', [0.75, 0.75, 0.75]);
%         plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
%         rectangle('Position', [lightOnsetTime 200 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
%         axis([xmin xmax ymin ymax])
%         xlabel('Time (s)');
%         ylabel(strcat("Normalized ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
%         title([obj.file ' (all) - Normalized'],'Interpreter','none');
%         hold off;
    
    
    % plotting figure with only Ch1 niceplot
        figure('name', strcat(obj.file,' (all) - niceplot - Normalized')); % naming figure file
        hold on;
        plot(x, dataPerSweepCh1,'Color', [0.75, 0.75, 0.75]);
        plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
        rectangle('Position', [lightOnsetTime 200 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
        axis([lightOnsetTime-0.01 lightOnsetTime+0.1 -2000 500]) % used to be -7000
        xlabel('Time (s)');
        ylabel(strcat("Normalized ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
        title([obj.file ' (all) - Normalized'],'Interpreter','none');
        hold off;
%     set(gca,'Visible','off')
    

 % save csv file with data 
        filename = strcat(obj.file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - normalized light evoked currents");
        fulldirectory = strcat(savefileto,'\',filename,'.csv');        
        dataInCellFormat = {};
        dataInCellFormat = num2cell(data);
        labeledData = cell2table(dataInCellFormat,'VariableNames',...
            {'sweep', 'seriesResistance', 'lightEvokedCurrent'});
        writetable(labeledData,fulldirectory);
        disp('I saved it')
        disp('Change directory if you want this saved elsewhere!')     
 
        
 % plot rs
        figure('name', strcat(obj.file,' (all) - rs')); % naming figure file
        plot(allSweeps, allRs,'-o');
        % plot lines marking 30% increase and 30% decrese in Rs compared to first
        % test pulse
        line([allSweeps(1) allSweeps(end)],[allRs(1)*0.7, allRs(1)*0.7],'Color','black','LineStyle','--')
        line([allSweeps(1) allSweeps(end)],[allRs(1)*1.3, allRs(1)*1.3],'Color','black','LineStyle','--')
        axis([allSweeps(1) allSweeps(end) 0 60])
        ylabel('Rs (M\Omega)');
        xlabel('Sweeps');
        title([obj.file ' rs'],'Interpreter','none');
        movegui('northeast');

end