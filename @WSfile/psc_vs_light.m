%{ 
DOCUMENTATION
Created: 2021 08 18
Author: PA

This function is used to analyze the amplitude of post-synaptic currents
(PSCs), either IPSCs, EPSCs, or ChR2-mediated currents.

The goal in mind when making this function is to analyse the release
probability of the connections between SNr and SNc.

This code is a mashup of firing_vs_light and normmono.

INPUTS explained:
    TO DO

INPUTS defaults:
    % Affects data analysis:
    lightStimCh = 2;
    discardedSweeps = [];
    discardedSweepsFromEnd = 1;
    inwardORoutward = 1;    % 1 (positive) is outward; -1 (negative) in inward
    baselineDurationInSeconds = 0.5;
    lightPulseAnalysisWindowInSeconds = 0.02;
    thresholdInDataPts = 10;
    rsTestPulseOnsetTime = 1;

    % Affects data display:
    ymax = 150;

    % Affects data saving:
    savefileto = 'D:\CORONAVIRUS DATA\From MATLAB';
    
OUTPUTS:
    TO DO

ASSUMPTIONS: 
    TO DO

TO DO:
    TO DO
%}

function psc_vs_light(obj)
%%%  USER INPUT ==================================================

% Affects data analysis:
lightStimCh = 2;
discardedSweeps = [];
discardedSweepsFromEnd = 0;
inwardORoutward = 1;    % 1 (positive) is outward; -1 (negative) in inward
baselineDurationInSeconds = 0.5;
lightPulseAnalysisWindowInSeconds = 0.02;
thresholdInDataPts = 10;
rsTestPulseOnsetTime = 1;

% Affects data display:
ymax = 150;

% Affects data saving:
savefileto = 'D:\CORONAVIRUS DATA\From MATLAB';


%% PREP - get info from file and create arrays ==================

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% getting info from file
samplingFrequency = obj.header.Acquisition.SampleRate;
fileName = obj.file;

% add 'subset' to fileName in case discardedSweeps is not empty
% to prevent overwritting of saved files with the whole dataset
if ~isempty(discardedSweeps)
    fileName = strcat(obj.file, '_subset');
end

% getting sweep numbers from file name
[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj); 

% checking for incomplete sweeps and not analyzing incomplete sweeps - to
% avoid this error: "Index exceeds the number of array elements (0)".   
if numel(fieldnames(obj.sweeps)) <= obj.header.NSweepsPerRun  
    lastSweepNumber = firstSweepNumber + numel(fieldnames(obj.sweeps)) -1 -discardedSweepsFromEnd;
    allSweeps = firstSweepNumber:lastSweepNumber;
end

% removing sweeps that should be discarded based on user input
for i=discardedSweeps
    allSweeps = allSweeps(allSweeps~=i);
end

% calculating variables based on user input
baselineDurationInDataPts = baselineDurationInSeconds * samplingFrequency;
lightPulseAnalysisWindowInDataPts = lightPulseAnalysisWindowInSeconds * samplingFrequency;
rsBaselineDataPointInterval = ((rsTestPulseOnsetTime-0.05)*samplingFrequency):(rsTestPulseOnsetTime*samplingFrequency);
rsFirstTransientDataPointInterval = (rsTestPulseOnsetTime*samplingFrequency):(rsTestPulseOnsetTime+0.0025)*samplingFrequency;
mouseNumber = getMouseNumber(obj);
experimentDate = getExperimentDate(obj);

% creating matrixes
lightEvokedCurrents = [];
allLightEvokedResponseLatencyInMilliSeconds = [];
allRs = [];
yBaselineSubAll = [];
lightEvokedCurrentsAllSweeps = [];
allLightEvokedResponseLatencyInMilliSecondsAllSweeps = [];
baselineCurrentAll = [];
data = [];



%% SWEEP BY SWEEP ANALYSIS ======================================

% get data from all sweeps in file
for sweepNumber = allSweeps
    
    % get light stim data
    [xch2,ych2] = obj.xy(sweepNumber, lightStimCh);      
    
    % get light stim parameters
    lightPulseStart = find(diff(ych2>1)>0);                                      % list of light pulse onset data points
    lightPulseEnd = find(diff(ych2<1)>0);                                        % list of light pulse offset data points
    lightOnsetTime = lightPulseStart(1)/samplingFrequency;                       % onset of first light pulse in seconds
    stimDur = (lightPulseEnd(end)-lightPulseStart(end))/samplingFrequency;       % duration of each light pulse in the train (s)
    stimInterval = (lightPulseStart(2)-lightPulseStart(1))/samplingFrequency;    % interval between each pulse (s)
    stimFreq = 1/stimInterval;                                                   % frequency of the light stim (Hz)
    lightDur = (lightPulseStart(end)-lightPulseStart(1))/samplingFrequency + stimInterval;  % duration of the whole light train stim (s)  
    
    % get raw data
    [x,y] = obj.xy(sweepNumber, 1);
    
    % baseline subtraction
    baselineStart = lightPulseStart(1) - baselineDurationInDataPts;
    yBaselineSub = y-mean(y(baselineStart:lightPulseStart(1)));
    baselineCurrent = mean(y(baselineStart:lightPulseStart(1)));
    
    % saving data for niceplot
    % y data for each sweep is in a column
    yBaselineSubAll = [yBaselineSubAll, yBaselineSub]; 
    baselineCurrentAll = [baselineCurrentAll, baselineCurrent];
    
%     % filter data
%     yFiltered = bandpass(y,[highpassThreshold lowpassThreshold],samplingFrequency);
%     
%     % find peaks or valleys based on user input
%     if peaksOrValleys == 'peaks'
%         [pks,locs,w,p] = findpeaks(yFiltered,x,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);
%     else
%         [pks,locs,w,p] = findpeaks(-yFiltered,x,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);
%     end

    % calculating series resistance
    rsBaselineCurrent = mean(y(rsBaselineDataPointInterval));
    rsTransientCurrent = min(y(rsFirstTransientDataPointInterval));
    dCurrent = rsTransientCurrent-rsBaselineCurrent;
    dVoltage = -5;
    seriesResistance = 1000*dVoltage/dCurrent; %mV/pA equals Gohm
    
    % put series resistance values from each sweep into a different column
    allRs = [allRs, seriesResistance];
    
    % clean up matrices for that will be used in the next loop
    lightEvokedCurrents = [];
    allLightEvokedResponseLatencyInMilliSeconds = [];

    % loop through all light pulses in the train
    for pulseOnset = lightPulseStart.'
        
        afterLightDataPoint = pulseOnset + lightPulseAnalysisWindowInDataPts;
            
        % get amplitude of lightEvokedCurrents
        % outward current is +1
        % inward current is 0 or -1
        if inwardORoutward == 1 
            lightEvokedCurrent = max(yBaselineSub(pulseOnset:afterLightDataPoint));
        else
            lightEvokedCurrent = min(yBaselineSub(pulseOnset:afterLightDataPoint));
        end
        
        % put all currents from a particular sweep in a row (each column is
        % a pulse)
        lightEvokedCurrents = [lightEvokedCurrents, lightEvokedCurrent];
        
        % get latency of lightEvokedCurrents
        %%% Rationale for finding light evoked response latency: I am looking
        % for a monotonic change of the signal for at least "x" data points 
        % (x=threshold), and I am selecting the first occurence of this monotonic change.       
        %%% The expected direction of the monotonic change is determined
        % by the optarg "inwardORoutward". If the value is 1, the function
        % will look for a monotonic increase (outward current), and if the
        % value is -1, the function will look for a monotonic decrease
        % (inward current).        
        %%% y(lightOnsetDataPoint:lightOffDataPoint) is the signal during light pulse        
        %%% zeros(threshold,1)+(inwardORoutward) is a vector with "threshold" 
        % number of rows containing -1 or 1, depending on the value of the
        % optarg "inwardORoutward"        
        %%% function diff calculates difference between data points        
        %%% function sign only keeps info about decrease (-1) or increase (1) in signal        
        %%% function conv2 performs convolution, looking for a sequence of
        % "threshold" points of value equal to the value stored in
        % "inwardORoutward" (-1 or 1). I think conv2 collapses all found 
        % elements into a single element of value "threshold" and moves 
        % forward in its search. So if threshold=20, inwardORoutward=-1,
        % and there are 40 elements of value -1 in sequence (one after the
        % other), conv2 will spit out a vector with 2 elements, both of
        % value 20.       
        %%% function find looks for the index of the elements equal to threshold in the output of the convolution        
        %%% function min looks for the min index (aka the first data point
        % that marks the start of a monotonic change in signal for "x"
        % points (x=threshold).
        lightEvokedResponseLatencyInDataPoints = min(find(conv2(sign(diff(y(pulseOnset:afterLightDataPoint))), zeros(thresholdInDataPts,1)+(inwardORoutward), 'valid')==thresholdInDataPts));
        
        % 1000 multiplication is conversion from seconds to milliseconds
        lightEvokedResponseLatencyInMilliSeconds = 1000*lightEvokedResponseLatencyInDataPoints/samplingFrequency;
        
        if isempty(lightEvokedResponseLatencyInMilliSeconds)
            lightEvokedResponseLatencyInMilliSeconds = NaN;
        end
        
        % put all latencies from a particular sweep in a row
        allLightEvokedResponseLatencyInMilliSeconds = [allLightEvokedResponseLatencyInMilliSeconds, lightEvokedResponseLatencyInMilliSeconds];        
        
    end
    
    % store all currents from a cell (each column is a pulse; each row is a sweep)
    lightEvokedCurrentsAllSweeps = [lightEvokedCurrentsAllSweeps; lightEvokedCurrents];
    allLightEvokedResponseLatencyInMilliSecondsAllSweeps = [allLightEvokedResponseLatencyInMilliSecondsAllSweeps; allLightEvokedResponseLatencyInMilliSeconds];
    
    % Data that will be exported (each sweep is a row)       
    data = [data; ...
        mouseNumber, ...
        experimentDate, ...
        sweepNumber, ...
        discardedSweepsFromEnd, ...
        inwardORoutward, ...
        baselineDurationInSeconds, ...
        lightPulseAnalysisWindowInSeconds, ...
        thresholdInDataPts, ...
        rsTestPulseOnsetTime, ...
        stimDur, ... 
        stimFreq, ...
        lightDur, ...
        seriesResistance, ...
        baselineCurrent];     
        
end

data = [data, lightEvokedCurrentsAllSweeps, allLightEvokedResponseLatencyInMilliSecondsAllSweeps];
    
xmin = lightOnsetTime-baselineDurationInSeconds;
xmax = lightOnsetTime+lightDur+baselineDurationInSeconds;


%% CELL Analysis

dataCell = [mouseNumber, ...
    experimentDate, ...
    firstSweepNumber, ...
    lastSweepNumber, ...
    length(allSweeps), ...
    mean(allRs), ...
    min(allRs), ...
    max(allRs), ...
    mean(baselineCurrentAll), ...
    min(baselineCurrentAll), ...
    max(baselineCurrentAll), ...
    mean(lightEvokedCurrentsAllSweeps), ...
    std(lightEvokedCurrentsAllSweeps), ...
    mean(allLightEvokedResponseLatencyInMilliSecondsAllSweeps, 'omitnan'), ...
    std(allLightEvokedResponseLatencyInMilliSecondsAllSweeps, 'omitnan')];


%% Color-blind diverging color scheme

brownToPurpleHex = {'#7f3b08','#b35806','#e08214','#fdb863','#fee0b6','#d8daeb','#b2abd2','#8073ac','#542788','#2d004b'};

color1rgb = [127, 59, 8];
color2rgb = [179, 88, 6];
color3rgb = [224, 130, 20];
color4rgb = [253, 184, 99];
color5rgb = [254, 224, 182];
color6rgb = [216, 218, 235];
color7rgb = [178, 171, 210];
color8rgb = [128, 115, 172];
color9rgb = [84, 39, 136];
color10rgb = [45, 1, 75];

brownToPurpleRgb = [color1rgb/255
                    color2rgb/255
                    color3rgb/255
                    color4rgb/255
                    color5rgb/255
                    color6rgb/255
                    color7rgb/255
                    color8rgb/255
                    color9rgb/255
                    color10rgb/255];

brownToPurpleRgbTransparent = [color1rgb/255, 0.25
                                color2rgb/255, 0.25
                                color3rgb/255, 0.25
                                color4rgb/255, 0.25
                                color5rgb/255, 0.25
                                color6rgb/255, 0.25
                                color7rgb/255, 0.25
                                color8rgb/255, 0.25
                                color9rgb/255, 0.25
                                color10rgb/255, 0.25];
                       

%% PLOT - light-evoked currents amplitude ==========================

f1=figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light - oPSCs amplitude all')); % naming figure file
hold on;
colororder(f1,brownToPurpleHex);
plot(lightEvokedCurrentsAllSweeps.');
colormap(brownToPurpleRgb(1:length(allSweeps),:))
c = colorbar('Ticks',[1,length(allSweeps)]);
caxis([1 length(allSweeps)]);
c.Label.String = 'Sweep #';
set(c, 'YDir', 'reverse');
line([0 60],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
axis([-inf inf -ymax ymax]);
title([fileName ' - psc_vs_light - oPSC amplitude all'],'Interpreter','none');
ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
xlabel('Light pulse');
hold off;
movegui('northwest');


%% PLOT - light-evoked currents latency ===========================

f2=figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light - oPSCs latency all')); % naming figure file
hold on;
colororder(f2,brownToPurpleHex);
plot(allLightEvokedResponseLatencyInMilliSecondsAllSweeps.','-o');
colormap(brownToPurpleRgb(1:length(allSweeps),:))
c = colorbar('Ticks',[1,length(allSweeps)]);
caxis([1 length(allSweeps)]);
c.Label.String = 'Sweep #';
set(c, 'YDir', 'reverse');
line([0 60],[5, 5],'Color','black','LineStyle','--')
line([0 60],[1, 1],'Color','red','LineStyle','--')
axis([-inf inf 0 10]);
title([fileName ' - psc_vs_light - oPSC latency all'],'Interpreter','none');
ylabel('Response Latency (ms)');
xlabel('Light pulse');
hold off;
movegui('north');


%% PLOT - rs =====================================================

figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light - Rs all')); % naming figure file
plot(allSweeps, allRs,'-o');
% plot lines marking 30% increase and 30% decrese in Rs compared to first
% test pulse
line([allSweeps(1) allSweeps(end)],[allRs(1)*0.7, allRs(1)*0.7],'Color','black','LineStyle','--')
line([allSweeps(1) allSweeps(end)],[allRs(1)*1.3, allRs(1)*1.3],'Color','black','LineStyle','--')
axis([allSweeps(1) inf 0 60])
ylabel('Rs (M\Omega)');
xlabel('Sweeps');
title([obj.file ' rs'],'Interpreter','none');
movegui('northeast');


%% PLOT - niceplot of all sweeps - COLORS ================================

% plotting niceplot     
figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light - Niceplot all colors'));
hold on;
for sweep = 1:size(yBaselineSubAll,2)
    plot(x, yBaselineSubAll(:,sweep),'Color',brownToPurpleRgbTransparent(sweep,:));
end
colormap(brownToPurpleRgb(1:sweep,:))
c = colorbar('Ticks',[1,sweep]);
caxis([1 sweep]);
c.Label.String = 'Sweep #';
set(c, 'YDir', 'reverse');
% plot(x, mean(yBaselineSubAll,2),'Color','black','LineWidth',1.5); 
line([xmin xmax],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
axis([xmin xmax -ymax ymax]);
xlabel('Time (s)');
ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
title([fileName ' - psc_vs_light - niceplot all'],'Interpreter','none');
set(gca,'Visible','off');
set(gcf,'Position',[1400 550 500 400]);

% adding light stim - individual pulses   
% note that this code will use light stim parameters from the last sweep!
% if light stim is not the same accross all sweeps, this will be
% misleading!
for nStim=1:length(lightPulseStart)
%     line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-ymax+100,-ymax+100],'Color',[0 0.4470 0.7410],'LineWidth',10)
    line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-inwardORoutward*(ymax/5),-inwardORoutward*(ymax/5)],'Color',[0 0.4470 0.7410],'LineWidth',10)
end

% adding scale bar
ymin = -ymax;
line([xmax-1 xmax],[ymax ymax],'Color','k')
line([xmax xmax],[ymax ymax-((ymax-ymin)/10)],'Color','k')
text(xmax-1, ymax-((ymax-ymin)/20), "1 s")
text(xmax-1, ymax-((ymax-ymin)/10), strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
hold off;
movegui('southeast');


%% PLOT - niceplot of all sweeps and/or MEAN - BLACK ================================

% plotting niceplot     
figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light - Niceplot black'));
hold on;
% plot(x, yBaselineSubAll,'Color',[0, 0, 0, 0.25]);
plot(x, mean(yBaselineSubAll,2),'Color','black','LineWidth',0.7); 
line([xmin xmax],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
axis([xmin xmax -ymax ymax]);
xlabel('Time (s)');
ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
title([fileName ' - psc_vs_light - niceplot all'],'Interpreter','none');
set(gca,'Visible','off');
set(gcf,'Position',[1400 550 500 400]);

% adding light stim - individual pulses   
% note that this code will use light stim parameters from the last sweep!
% if light stim is not the same accross all sweeps, this will be
% misleading!
for nStim=1:length(lightPulseStart)
    line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-inwardORoutward*(ymax/5),-inwardORoutward*(ymax/5)],'Color',[0 0.4470 0.7410],'LineWidth',10)
end

% adding light stim - train
% note that this code will use light stim parameters from the last sweep!
% if light stim is not the same accross all sweeps, this will be
% misleading!
% line([lightOnsetTime,lightOnsetTime+lightDur],[-ymax+50,-ymax+50],'Color',[0 0.4470 0.7410],'LineWidth',10)

% adding scale bar
ymin = -ymax;
line([xmax-1 xmax],[ymax ymax],'Color','k')
line([xmax xmax],[ymax ymax-((ymax-ymin)/10)],'Color','k')
text(xmax-1, ymax-((ymax-ymin)/20), "1 s")
text(xmax-1, ymax-((ymax-ymin)/10), strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
hold off;
movegui('south');


%% PLOT - subtracted baseline current =====================================================

figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light - baseline current')); % naming figure file
plot(allSweeps, baselineCurrentAll,'-o');
axis([allSweeps(1) inf -300 300])
ylabel('Subtracted Baseline Current (pA)');
xlabel('Sweeps');
title([obj.file ' rs'],'Interpreter','none');
movegui('southwest');


%% PLOT - oPSC amplitude mean +- SD ================================

figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light - mean and SD'));
hold on;
errorbar(mean(lightEvokedCurrentsAllSweeps), std(lightEvokedCurrentsAllSweeps), 'color', [0 0 0], 'CapSize', 0);
% errorbar(mean(lightEvokedCurrentsAllSweeps), std(lightEvokedCurrentsAllSweeps),'-o', 'MarkerFaceColor', 'k', 'color', [0 0 0], 'CapSize', 0);
line([0 60],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
axis([-inf inf -ymax ymax]);
ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
title([fileName ' - psc_vs_light - mean and SD'],'Interpreter','none');
xlabel('Light pulse');
hold off;
movegui('northwest');


%% EXPORTING XLS files ==========================================

% create cell array with strings for naming the amplitude and latency of the oPSC for each light pulse
oPSCvariableNamesAmplitude = cell(1,length(lightEvokedCurrents));
oPSCvariableNamesLatency = cell(1,length(allLightEvokedResponseLatencyInMilliSeconds));
for lightPulse = 1:length(lightEvokedCurrents)
    oPSCvariableNamesAmplitude(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(pA)')};
    oPSCvariableNamesLatency(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(ms)')};
end

% stores sweep by sweep data
filename = strcat(fileName, '_', analysisDate, " - psc_vs_light - sweep_by_sweep");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(data);
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...      
    {'mouse', ...
    'date', ...
    'sweep', ...
    'discardedSweepsFromEnd', ...
    'inward(-1)ORoutward(1)', ...
    'baselineDurationInSeconds', ...
    'lightPulseAnalysisWindowInSeconds', ...
    'thresholdInDataPts', ...        
    'rsTestPulseOnsetTime', ...    
    'lightPulseDur(s)', ... 
    'lightStimFreq(Hz)', ...
    'lightDur(s)', ...
    'seriesResistance(Mohm)', ...
    'baselineCurrent(pA)', ...
    oPSCvariableNamesAmplitude{:}, ...
    oPSCvariableNamesLatency{:}});    
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the sweep_by_sweep xls file')


% create cell array with strings for naming the mean and std of the oPSC for each light pulse
oPSCvariableNamesAmplitudeAVG = cell(1,length(lightEvokedCurrents));
oPSCvariableNamesAmplitudeSTD = cell(1,length(lightEvokedCurrents));
oPSCvariableNamesLatencyAVG = cell(1,length(allLightEvokedResponseLatencyInMilliSeconds));
oPSCvariableNamesLatencySTD = cell(1,length(allLightEvokedResponseLatencyInMilliSeconds));
for lightPulse = 1:length(lightEvokedCurrents)
    oPSCvariableNamesAmplitudeAVG(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(pA)AVG')};
    oPSCvariableNamesAmplitudeSTD(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(pA)STD')};
    oPSCvariableNamesLatencyAVG(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(ms)AVG')};
    oPSCvariableNamesLatencySTD(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(ms)STD')};
end

% stores cell data
filename = strcat(fileName, '_', analysisDate, " - psc_vs_light - cell");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(dataCell);
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...      
    {'mouse', ...
    'date', ...
    'firstSweep', ...
    'lastSweep', ...
    'nSweeps', ...     
    'seriesResistanceAGV(Mohm)', ...
    'seriesResistanceMIN(Mohm)', ...
    'seriesResistanceMAX(Mohm)', ...   
    'baselineCurrentAVG(pA)', ...
    'baselineCurrentMIN(pA)', ...
    'baselineCurrentMAX(pA)', ...    
    oPSCvariableNamesAmplitudeAVG{:}, ...
    oPSCvariableNamesAmplitudeSTD{:}, ...
    oPSCvariableNamesLatencyAVG{:}, ...
    oPSCvariableNamesLatencySTD{:}});    
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the cell xls file')


end    