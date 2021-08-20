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
    TO DO
    
OUTPUTS:
    TO DO

ASSUMPTIONS: 
    TO DO

TO DO:
    TO DO
%}

function psc_vs_light(obj)
%%%  USER INPUT ==================================================

discardedSweeps = [];
discardedSweepsFromEnd = 1;
inwardORoutward = 1;
baselineDurationInSeconds = 0.5;
lightPulseAnalysisWindowInSeconds = 0.02;
thresholdInDataPts = 10;
rsTestPulseOnsetTime = 1;
ymax = 150;


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



%% SWEEP BY SWEEP ANALYSIS ======================================

% get data from all sweeps in file
for sweepNumber = allSweeps
    
    % get light stim data
    [xch2,ych2] = obj.xy(sweepNumber, 2);      
    
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
    
    % saving data for niceplot
    % y data for each sweep is in a column
    yBaselineSubAll = [yBaselineSubAll, yBaselineSub]; 
    
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
    
end
    
xmin = lightOnsetTime-baselineDurationInSeconds;
xmax = lightOnsetTime+lightDur+baselineDurationInSeconds;


%% PLOT - light-evoked currents amplitude ==========================

figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light - oPSCs amplitude all')); % naming figure file
hold on;
plot(lightEvokedCurrentsAllSweeps.','-o');
line([0 60],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
axis([-inf inf -ymax ymax]);
title([fileName ' - psc_vs_light - oPSC amplitude all'],'Interpreter','none');
ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
xlabel('Light pulse');
hold off;


%% PLOT - light-evoked currents latency ===========================

figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light - oPSCs latency all')); % naming figure file
hold on;
plot(allLightEvokedResponseLatencyInMilliSecondsAllSweeps.','-o');
line([0 60],[5, 5],'Color','black','LineStyle','--')
line([0 60],[1, 1],'Color','red','LineStyle','--')
axis([-inf inf 0 10]);
title([fileName ' - psc_vs_light - oPSC latency all'],'Interpreter','none');
ylabel('Response Latency (ms)');
xlabel('Light pulse');
hold off;


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


%% PLOT - niceplot of all sweeps ================================

% plotting niceplot     
figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light - Niceplot all'));
hold on;
line([xmin xmax],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
plot(x, yBaselineSubAll,'Color',[0, 0, 0, 0.25]);
plot(x, mean(yBaselineSubAll,2),'Color','black','LineWidth',1.5); 
axis([xmin xmax -ymax ymax]);
xlabel('Time (s)');
ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
title([fileName ' - psc_vs_light - niceplot all'],'Interpreter','none');
% set(gca,'Visible','off');
set(gcf,'Position',[1400 550 500 400]);

% adding light stim - individual pulses   
% note that this code will use light stim parameters from the last sweep!
% if light stim is not the same accross all sweeps, this will be
% misleading!
for nStim=1:length(lightPulseStart)
    line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-ymax+100,-ymax+100],'Color',[0 0.4470 0.7410],'LineWidth',10)
end

% adding light stim - train
% note that this code will use light stim parameters from the last sweep!
% if light stim is not the same accross all sweeps, this will be
% misleading!
line([lightOnsetTime,lightOnsetTime+lightDur],[-ymax+50,-ymax+50],'Color',[0 0.4470 0.7410],'LineWidth',10)

% adding scale bar
ymin = -ymax;
line([xmax-1 xmax],[ymax ymax],'Color','k')
line([xmax xmax],[ymax ymax-((ymax-ymin)/10)],'Color','k')
text(xmax-1, ymax-((ymax-ymin)/20), "1 s")
text(xmax-1, ymax-((ymax-ymin)/10), strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
hold off;

        
end    