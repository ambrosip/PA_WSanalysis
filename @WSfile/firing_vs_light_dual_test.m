%{ 
DOCUMENTATION
Created: 2021 09 13
Author: PA

This function is used to test the input parameters to be used in the 
function "firing_vs_light_dual". Use it with caution - it automatically
identifies the light stim parameters but the data collected is manually
adjusted to fit a 9s analysis window (3 s before/during/after light stim).

It extends the function firing_vs_light_test.

INPUTS explained:
    - discardedSweeps: specific sweeps that should NOT be analyzed due to
    artifacts. If I want to discard sweep 0024, just write 24.

    - discardedSweepsFromEnd: number of sweeps from the end that will not be
    analyzed. Used to remove incomplete sweeps or sweeps with artifacts. 

    - peaksOrValleys: find APs based on found peaks ('peaks') or valleys
    ('v').

    - highpassThreshold: bandpass filter parameter - attenuates low
    frequency noise below this frequency. Removes baseline drift.

    - lowpassThreshold: bandpass filter parameter - attenuates high
    frequency noise above this frequency. Removes fast artifacts (such as
    light stim onset/offset and outflow bubbles). Should be > 1000 to avoid
    sinusoidal artifacts. Use 2000 for striatum cells.

    - minPeakHeight: amplitude threshold for finding APs - detect peaks or
    valleys above/below this amplitude threshold.

    - minPeakDistance: min distance in seconds between found peaks. Used to
    ignore fast artifacts that cross the minPeakHeight threshold. Note that
    this input will cap the max firing frequency output. In other words: to
    further discard high frequency noise, only look for peaks that are this
    time interval apart (note that this interval caps the max signal
    frequency you can detect).

    - ymax: for illustration purposes, use this value as the max range for
    the y-axis when plotting current vs time.

    - ymaxhist: for illustration purposes, use this value as the max range
    for the y-axis in the firing rate histogram (max instantaneous firing
    frequency).

    - zoomWindow: time before/after LightOnsetTime that will be shown in
    the zoomed plot.

    - ymaxIsiCV: for illustration purposes, use this value as the max range
    for the y-axis in the ISI histogram.

    - savefileto: save CSV files to this folder.

INPUTS defaults:
    discardedSweeps = [];
    discardedSweepsFromEnd = 0;
    peaksOrValleys = 'v';   
    highpassThreshold = 100;
    lowpassThreshold = 1500;    
    minPeakHeight = 15;         
    minPeakDistance = 0.025;    

    preAPinSeconds = 0.005;            
    postAPinSeconds = 0.01;           
    preAPbaselineDurationSeconds = 0.002;
    ddyValleyThreshold = 30;
    ddyPeakThreshold = 10;

    ymax = 75;
    ymaxhist = 15;
    ZoomWindow = 0.25;
    ymaxIsiCV = 150;
    
OUTPUTS:
    Figs

ASSUMPTIONS: 
    - Recording was done in cell attached or loose cell mode in VClamp -
    bandpass filter is appropriate. 

TO DO:
    - complete documentation
%}

function firing_vs_light_dual_test(obj)


%%  USER INPUT ==================================================

% Affects data analysis - Finding APs:
discardedSweeps = [];
discardedSweepsFromEnd = 0;
peaksOrValleys = 'v';   
highpassThreshold = 100;
lowpassThreshold = 1500;    
minPeakHeight = 25;         
minPeakDistance = 0.005;    

% Affects data analysis - AP shape:
preAPinSeconds = 0.005;            
postAPinSeconds = 0.01;           
preAPbaselineDurationSeconds = 0.002;
ddyValleyThreshold = 50;
ddyPeakThreshold = 30;
  
% Affects data display: 
ymax = 100;
ymaxhist = 100;
zoomWindow = 3;
ymaxIsiCV = 150;


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

% creating matrixes/arrays that will be filled for AP shape
xSubset = [];
ySubsetForAPshape = [];
yFilteredSubset = [];
xSubsetAll = [];
ySubsetAll = [];
yFilteredSubsetAll = [];
nAPtotal = 0;


%% SWEEP BY SWEEP ANALYSIS ======================================

% get data from all sweeps in file
for sweepNumber = allSweeps
    
    % get raw data
    [x,y] = obj.xy(sweepNumber, 1);
    
    % filter data
    yFiltered = bandpass(y,[highpassThreshold lowpassThreshold],samplingFrequency);
    
    % find peaks or valleys based on user input
    if peaksOrValleys == 'peaks'
        [pks,locs,w,p] = findpeaks(yFiltered,x,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);
    else
        [pks,locs,w,p] = findpeaks(-yFiltered,x,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);
    end

    % get light stim data for channel 2
    [xch2,ych2] = obj.xy(sweepNumber, 2);      
    
    % get light stim parameters for channel 2
    lightPulseStartCh2 = find(diff(ych2>1)>0);
    lightPulseEndCh2 = find(diff(ych2<1)>0);
    lightOnsetTimeCh2 = lightPulseStartCh2(1)/samplingFrequency;                                        % in seconds
    stimDurCh2 = (lightPulseEndCh2(end)-lightPulseStartCh2(end))/samplingFrequency;                     % duration of each light pulse in the train (s)
    stimIntervalCh2 = (lightPulseStartCh2(2)-lightPulseStartCh2(1))/samplingFrequency;                  % interval between each pulse (s)
    stimFreqCh2 = 1/stimIntervalCh2;                                                                    % frequency of the light stim (Hz)
    lightDurCh2 = (lightPulseStartCh2(end)-lightPulseStartCh2(1))/samplingFrequency + stimIntervalCh2;  % duration of the whole light train stim (s)    
    
    
    % get light stim data for channel 3
    [xch3,ych3] = obj.xy(sweepNumber, 3);      
    
    % get light stim parameters for channel 3
    % I commented out some stuff cuz so far ch3 is a single, continuous pulse!
    lightPulseStartCh3 = find(diff(ych3>1)>0);
    lightPulseEndCh3 = find(diff(ych3<1)>0);
    lightOnsetTimeCh3 = lightPulseStartCh3(1)/samplingFrequency;                                        % in seconds
    stimDurCh3 = (lightPulseEndCh3(end)-lightPulseStartCh3(end))/samplingFrequency;                     % duration of each light pulse in the train (s)
%     stimIntervalCh3 = (lightPulseStartCh3(2)-lightPulseStartCh3(1))/samplingFrequency;                  % interval between each pulse (s)
%     stimFreqCh3 = 1/stimIntervalCh3;                                                                    % frequency of the light stim (Hz)
%     lightDurCh3 = (lightPulseStartCh3(end)-lightPulseStartCh3(1))/samplingFrequency + stimIntervalCh3;  % duration of the whole light train stim (s)  
    %----------------------------------------------------------------
    
    
    % Getting peak/valley amplitudes pre, during and post light
    % to do manual quality control of found peaks later
    if peaksOrValleys == 'peaks'
        pksPreLightBin1 = pks;
        pksPreLightBin2 = pks;
        pksPreLightBin3 = pks;
        pksDuringLightBin1 = pks;
        pksDuringLightBin2 = pks;
        pksDuringLightBin3 = pks;
        pksPostLightBin1 = pks;
        pksPostLightBin2 = pks;     
        pksPostLightBin3 = pks;
    else
        pksPreLightBin1 = -pks;
        pksPreLightBin2 = -pks;
        pksPreLightBin3 = -pks;
        pksDuringLightBin1 = -pks;
        pksDuringLightBin2 = -pks;
        pksDuringLightBin3 = -pks;
        pksPostLightBin1 = -pks;
        pksPostLightBin2 = -pks;     
        pksPostLightBin3 = -pks;
    end
    %----------------------------------------------------------------
    
    
    % Getting timestamps pre, during, and post light  
        % locsPreLight will be different from locsBaseline: locsPreLight is
        % the baseline immediately before the light stim (from
        % lightOnsetTime-lightDur to lightOnsetTime), while locsBaseline is
        % the full baseline (from 0s to lightOnsetTime).
    % first, store locs (timestamps) in new variables 
    locsPreLightBin1 = locs;
    locsPreLightBin2 = locs;
    locsPreLightBin3 = locs;
    locsDuringLightBin1 = locs;
    locsDuringLightBin2 = locs;
    locsDuringLightBin3 = locs;
    locsPostLightBin1 = locs;
    locsPostLightBin2 = locs;
    locsPostLightBin3 = locs;

    
    % let's go bin by bin, starting with preLight
    % ASSUMPTION: Light stim on Ch2 is longer than light stim on Ch3
    % locsPreLightBin1 is the first bin. It starts 3s before the onset of
    % light on Ch2 (and 4s before the onset of light on Ch3).
    indicesToDelete = find(locs<lightOnsetTimeCh2-3 | locs>=lightOnsetTimeCh2-2);
    locsPreLightBin1(indicesToDelete) = [];
    pksPreLightBin1(indicesToDelete) = [];
    
    indicesToDelete = find(locs<lightOnsetTimeCh2-2 | locs>=lightOnsetTimeCh2-1);
    locsPreLightBin2(indicesToDelete) = [];
    pksPreLightBin2(indicesToDelete) = [];
    
    indicesToDelete = find(locs<lightOnsetTimeCh2-1 | locs>=lightOnsetTimeCh2);
    locsPreLightBin3(indicesToDelete) = [];
    pksPreLightBin3(indicesToDelete) = [];
  
    % now duringLight
    indicesToDelete = find(locs<lightOnsetTimeCh2 | locs>(lightOnsetTimeCh2+1));
    locsDuringLightBin1(indicesToDelete) = [];
    pksDuringLightBin1(indicesToDelete) = [];
    
    indicesToDelete = find(locs<lightOnsetTimeCh2+1 | locs>(lightOnsetTimeCh2+2));
    locsDuringLightBin2(indicesToDelete) = [];
    pksDuringLightBin2(indicesToDelete) = [];
    
    indicesToDelete = find(locs<lightOnsetTimeCh2+2 | locs>(lightOnsetTimeCh2+3));
    locsDuringLightBin3(indicesToDelete) = [];
    pksDuringLightBin3(indicesToDelete) = [];
    
    % now postLight
    indicesToDelete = find(locs<(lightOnsetTimeCh2+lightDurCh2) | locs>(lightOnsetTimeCh2+lightDurCh2+1));
    locsPostLightBin1(indicesToDelete) = [];
    pksPostLightBin1(indicesToDelete) = [];
   
    indicesToDelete = find(locs<(lightOnsetTimeCh2+lightDurCh2+1) | locs>(lightOnsetTimeCh2+lightDurCh2+2));
    locsPostLightBin2(indicesToDelete) = [];
    pksPostLightBin2(indicesToDelete) = [];
    
    indicesToDelete = find(locs<(lightOnsetTimeCh2+lightDurCh2+2) | locs>(lightOnsetTimeCh2+lightDurCh2+3));
    locsPostLightBin3(indicesToDelete) = [];
    pksPostLightBin3(indicesToDelete) = [];
    %----------------------------------------------------------------
    
    
    % Collecting AP shape data  
    % getting all of the APs timestamps prior to light stim (aka full
    % baseline, from 0s to lightOnsetTime)
    locsBaseline = locs;
    indicesToDelete = find(locs >= lightOnsetTimeCh2);
    locsBaseline(indicesToDelete) = [];
    
    % Store total number of APs found in complete baseline period
    nAP = size(locsBaseline,1);
    
    % Store data surrounding each AP
    % Exclude first and last AP to avoid error - MATLAB will complain if
    % the required time window around each AP is beyond window of data
    % available.
    for iAP = 2:(nAP-1)
                
        % convert time points to data points
        xminForAPshape = round(samplingFrequency * (locsBaseline(iAP) - preAPinSeconds));
        xmaxForAPshape = round(samplingFrequency * (locsBaseline(iAP) + postAPinSeconds));
        xbaselineForAPshape = round(xminForAPshape + samplingFrequency*preAPbaselineDurationSeconds);
        
        % make data rows  
        ySubsetForAPshape = transpose(y(xminForAPshape:xmaxForAPshape));
        
        % calculate mean y during user-defined baseline
        ySubsetBaseline = mean(y(xminForAPshape:xbaselineForAPshape));      
        
        % calculate baseline subtracted y
        ySubsetBaselineSubtracted = ySubsetForAPshape - ySubsetBaseline;
        
        % Store all data
        ySubsetAll = [ySubsetAll; ySubsetBaselineSubtracted];
        
    end   
    
    nAPtotal = nAPtotal + nAP;
    %----------------------------------------------------------------
    
    
    % PLOT firing histogram for each sweep (1/ISI)
    sweepDuration = obj.header.Acquisition.Duration;
    sweepTime=0;
    inverseISIperSecBin=[];
    
    % naming figure file    
    figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_dual_test - ISI histogram'));      
    
        % Ch1 (recorded data) histogram
        subplot(3,1,1)         
        while sweepTime <= sweepTime+1 & sweepTime < sweepDuration
            indicesToDelete = find(locs<sweepTime | locs>=sweepTime+1);
            locsDuringSweepTime = locs;
            locsDuringSweepTime(indicesToDelete) = [];
            inverseISIperSecBin=[inverseISIperSecBin 1/mean(diff(locsDuringSweepTime))];
            sweepTime=sweepTime+1;
        end
        bar(0:length(inverseISIperSecBin)-1,inverseISIperSecBin,1);
        axis([-inf inf 0 ymaxhist]);            
        xlabel('Bins (1 s long)');
        ylabel('1/ISI (Hz)');
        title([' (' num2str(sweepNumber) ') - ISI hist'],'Interpreter','none');        

        % Ch2 (light stim) plot
        subplot(3,1,2)   
        plot(xch2,ych2);
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([0 30 yminhere ymaxhere])
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));        
        
        % Ch3 (light stim) plot
        subplot(3,1,3)   
        plot(xch3,ych3);
        yminhere = min(ych3)-5;
        ymaxhere = max(ych3)+5;
        axis([0 30 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(3), ' (', obj.header.Acquisition.AnalogChannelUnits(3), ')'));        
        set(gcf,'Position',[50 350 280 630]) % top left corner
    %----------------------------------------------------------------        
        
    
    % PLOT firing histogram for each sweep (#APs)
    sweepDuration = obj.header.Acquisition.Duration;
    sweepTime=0;
    nAPperSecBin=[];
    
    % naming figure file    
    figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_dual_test - nAPs histogram'));      
    
        % Ch1 (recorded data) histogram
        subplot(3,1,1)         
        while sweepTime <= sweepTime+1 & sweepTime < sweepDuration
            indicesToDelete = find(locs<sweepTime | locs>=sweepTime+1);
            locsDuringSweepTime = locs;
            locsDuringSweepTime(indicesToDelete) = [];
            nAPperSecBin = [nAPperSecBin length(locsDuringSweepTime)];
            sweepTime = sweepTime+1;
        end
        bar(0:length(nAPperSecBin)-1,nAPperSecBin,1);
        axis([-inf inf 0 ymaxhist]);            
        xlabel('Bins (1 s long)');
        ylabel('#APs (Hz)');
        title([' (' num2str(sweepNumber) ') - nAPs hist'],'Interpreter','none');        

        % Ch2 (light stim) plot
        subplot(3,1,2)   
        plot(xch2,ych2);
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([0 30 yminhere ymaxhere])
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));        
        
        % Ch3 (light stim) plot
        subplot(3,1,3)   
        plot(xch3,ych3);
        yminhere = min(ych3)-5;
        ymaxhere = max(ych3)+5;
        axis([0 30 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(3), ' (', obj.header.Acquisition.AnalogChannelUnits(3), ')'));        
        set(gcf,'Position',[350 350 280 630])
    %----------------------------------------------------------------  
        
    
    % PLOT found peaks for each sweep    
    hzPreLightBin3 = length(locsPreLightBin3);
    hzDuringLightBin1 = length(locsDuringLightBin1);
    hzDuringLightBin2 = length(locsDuringLightBin2);   
    
    figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_dual_test - found peaks')); 
        
        % Ch1 subplot (recorded data)
        subplot(3,1,1);
        plot(x,yFiltered)
        hold on;
        plot(locsPreLightBin3,pksPreLightBin3,'o','color','green');
        plot(locsDuringLightBin1,pksDuringLightBin1,'o','color','blue');
        plot(locsDuringLightBin2,pksDuringLightBin2,'o','color','red');
        hold off;
        axis([lightOnsetTimeCh2-4 lightOnsetTimeCh2+lightDurCh2+4 -inf ymax])
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([obj.file ' (' num2str(sweepNumber) ') - peaks'],'Interpreter','none');
        text(lightOnsetTimeCh2 + lightDurCh2/4,0.9*ymax,num2str(round(hzDuringLightBin1,2)),'color','blue');
        text(lightOnsetTimeCh2 - 3,0.9*ymax,num2str(round(hzPreLightBin3,2)),'color','green');
        text(lightOnsetTimeCh2 + lightDurCh2 + 2,0.9*ymax,num2str(round(hzDuringLightBin2,2)),'color','red');
        
        % Ch2 subplot (light stim command)
        subplot(3,1,2);
        plot(xch2,ych2)
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([lightOnsetTimeCh2-4 lightOnsetTimeCh2+lightDurCh2+4 yminhere ymaxhere])
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));
        
        % Ch3 subplot (light stim command)
        subplot(3,1,3);
        plot(xch3,ych3)
        yminhere = min(ych3)-5;
        ymaxhere = max(ych3)+5;
        axis([lightOnsetTimeCh2-4 lightOnsetTimeCh2+lightDurCh2+4 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(3), ' (', obj.header.Acquisition.AnalogChannelUnits(3), ')'));
        set(gcf,'Position',[650 350 500 600])
    %----------------------------------------------------------------         
    
    
%     % PLOT zoomed in niceplot (start of light stim)
%     figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_dual_test - zoom start')); 
%     
%         % main plot
%         plot(x,yFiltered,'k','LineWidth',1)
%         xmin = lightOnsetTimeCh2-zoomWindow;
%         xmax = lightOnsetTimeCh2+zoomWindow;
%         ymin = -ymax;
%         axis([xmin xmax ymin ymax]);
%         title([obj.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
%         set(gca,'Visible','off');
%         set(gcf,'Position',[50 50 500 400])
%         
%         % adding light stim Ch2
%         interStimInterval = 1/stimFreqCh2;
%         postStimInterval = interStimInterval - stimDurCh2;
%         for nStim = 1 : stimFreqCh2 * xmax
%             startStim = lightOnsetTimeCh2 + (nStim - 1) * interStimInterval;
%             line([startStim, startStim + stimDurCh2],[ymax,ymax],'Color','blue','LineWidth',10)
%         end
%         
%         % adding light stim Ch3
%         line([lightOnsetTimeCh3, lightOnsetTimeCh3 + stimDurCh3],[ymin,ymin],'Color','red','LineWidth',10)
% %         interStimInterval = 1/stimFreqCh3;
% %         postStimInterval = interStimInterval - stimDurCh3;
% %         for nStim = 1 : stimFreqCh3 * xmax
% %             startStim = lightOnsetTimeCh3 + (nStim - 1) * interStimInterval;
% %             line([startStim, startStim + stimDurCh3],[ymin,ymin],'Color','red','LineWidth',10)
% %         end
%         
%         % adding scale bar
%         line([xmin,xmin+zoomWindow],[ymin,ymin],'Color','k')
%         line([xmin,xmin],[ymin,ymin+((ymax-ymin)/10)],'Color','k')
%         text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/30),strcat(num2str(zoomWindow*1000)," ms"))
%         text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/12),strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
    %---------------------------------------------------------------- 
     
    
%     % PLOT zoomed in niceplot (end of light stim) 
%     figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_dual_test - zoom end')); 
%     
%         % main plot
%         plot(x,yFiltered,'k','LineWidth',1)
%         xmin=lightOnsetTimeCh2+lightDurCh2-zoomWindow;
%         xmax=lightOnsetTimeCh2+lightDurCh2+zoomWindow;
%         ymin=-ymax;
%         axis([xmin xmax ymin ymax]);
%         title([obj.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
%         set(gca,'Visible','off');
%         set(gcf,'Position',[600 50 500 400])
%         
%         % adding light stim
%         interStimInterval = 1/stimFreqCh2;
%         postStimInterval = interStimInterval-stimDurCh2;
%         for nStim=1:stimFreqCh2*zoomWindow
%             startStim = lightOnsetTimeCh2+lightDurCh2-zoomWindow + (nStim - 1) * interStimInterval;
%             line([startStim,startStim+stimDurCh2],[ymax,ymax],'Color',[0 0.4470 0.7410],'LineWidth',10)
%         end
%         
%         % adding scale bar
%         line([xmin,xmin+zoomWindow],[ymin,ymin],'Color','k')
%         line([xmin,xmin],[ymin,ymin+((ymax-ymin)/10)],'Color','k')
%         text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/30),strcat(num2str(zoomWindow*1000)," ms"))
%         text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/12),strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
    %----------------------------------------------------------------        
    
    
    % PLOT one figure per sweep that shows whole sweep (raw)    
    figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_dual_test - raw'));   
        
        subplot(3,1,1)
        plot(x,y);
        axis([-inf inf -inf inf]);           
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([' (' num2str(sweepNumber) ') - raw'],'Interpreter','none');
      
        subplot(3,1,2) 
        plot(xch2,ych2);
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([0 30 yminhere ymaxhere])
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));        
        
        subplot(3,1,3) 
        plot(xch3,ych3);
        yminhere = min(ych3)-5;
        ymaxhere = max(ych3)+5;
        axis([0 30 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(3), ' (', obj.header.Acquisition.AnalogChannelUnits(3), ')'));        
        set(gcf,'Position',[1200 350 280 630])
    %----------------------------------------------------------------    

    
    % PLOT one figure per sweep that shows whole sweep (filtered)    
    figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_dual_test - filtered'));  
        
        subplot(3,1,1)
        plot(x,yFiltered);
        axis([-inf inf -inf ymax]);           
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([' (' num2str(sweepNumber) ') - filtered'],'Interpreter','none');
      
        subplot(3,1,2) 
        plot(xch2,ych2);
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([0 30 yminhere ymaxhere])
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));        
        
        subplot(3,1,3) 
        plot(xch3,ych3);
        yminhere = min(ych3)-5;
        ymaxhere = max(ych3)+5;
        axis([0 30 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(3), ' (', obj.header.Acquisition.AnalogChannelUnits(3), ')'));        
        set(gcf,'Position',[1200 50 280 630])
    %----------------------------------------------------------------              
end


%% CELL ANALYSIS - AP shape =====================================

% calculate average AP shape/trace
avgAP = mean(ySubsetAll);
%----------------------------------------------------------------

% create x axis for plotting AP shape
xSubset = 1000*linspace(0,(preAPinSeconds + postAPinSeconds), length(ySubsetForAPshape));
%----------------------------------------------------------------

% find AP peak and valley
avgAPpeakInDataPoints = find(avgAP==max(avgAP));
avgAPvalleyInDataPoints = find(avgAP==min(avgAP));
avgAPpeakInMilliSeconds = xSubset(avgAPpeakInDataPoints);
avgAPvalleyInMilliSeconds = xSubset(avgAPvalleyInDataPoints);
%----------------------------------------------------------------

% calculating derivatives and creating xAxis to plot derivatives. 
% for each derivative, the xAxis length decreases by 1 point.
xForDy = xSubset;
xForDy(end) = [];   % remove last point from xAxis
dy = diff(avgAP)./diff(xSubset);
ddy = diff(dy)./diff(xForDy);
xForDdy = xForDy;
xForDdy(end) = [];  % remove last point from xAxis
%----------------------------------------------------------------

% Find AP ONset based on 2nd derivative - ddy == min of first valley
[pks1,locs1,w,p] = findpeaks(-ddy,xForDdy,'MinPeakHeight',ddyValleyThreshold);
ddyBasedOnsetInMilliSeconds = locs1(1);       
ddyBasedOnsetInDataPoints = round(ddyBasedOnsetInMilliSeconds*(samplingFrequency/1000));
%----------------------------------------------------------------

% Find AP OFFset based on 2nd derivative - ddy == zero after last peak
[pks2,locs2,w,p] = findpeaks(ddy,xForDdy,'MinPeakHeight',ddyPeakThreshold);
ddyLastPeakInMilliSeconds = locs2(end);     
ddyLastPeakInDataPoints = round(ddyLastPeakInMilliSeconds*(samplingFrequency/1000));
ddyLastValleyInMilliSeconds = locs1(end);     
ddyLastValleyInDataPoints = round(ddyLastValleyInMilliSeconds*(samplingFrequency/1000));
ddyAfterLastPeakOrValley = ddy;
ddyBasedOffsetInDataPoints = [];

% if last peak is before last valley, look for when ddy crosses 0 after last valley
if ddyLastPeakInMilliSeconds < ddyLastValleyInMilliSeconds
    ddyAfterLastPeakOrValley(1:ddyLastValleyInDataPoints) = [];
    ddyAfterLastPeakOrValleyInDataPoints = ddyLastValleyInDataPoints;
    ddyCrossesZeroPt = find(diff(sign(ddyAfterLastPeakOrValley))>0, 1);
    ddyBasedOffsetInDataPoints = ddyLastValleyInDataPoints + ddyCrossesZeroPt + 1;
% if last peak is after last valley, look for when ddy crosses 0 after last peak
else
    ddyAfterLastPeakOrValley(1:ddyLastPeakInDataPoints) = [];
    ddyAfterLastPeakOrValleyInDataPoints = ddyLastPeakInDataPoints;
    ddyCrossesZeroPt = find(diff(sign(ddyAfterLastPeakOrValley))<0, 1);
    ddyBasedOffsetInDataPoints = ddyLastPeakInDataPoints + ddyCrossesZeroPt + 1;
end 

ddyBasedOffsetInMilliSeconds = xForDdy(ddyBasedOffsetInDataPoints);
ddyAfterLastPeakOrValleyInMilliSeconds = xForDdy(ddyAfterLastPeakOrValleyInDataPoints);
%----------------------------------------------------------------

% Find AP offset based on avgAP==0 after peak
% Assumes that AP valley precedes the AP peak
avgAPafterPeak = avgAP;
avgAPafterPeak(1:avgAPpeakInDataPoints) = [];
avgAPoffsetInDataPoints = avgAPpeakInDataPoints + find(round(avgAPafterPeak)==0, 1);
avgAPoffsetInMilliSeconds = xSubset(avgAPoffsetInDataPoints);

% To avoid errors in the CSV in case matlab fails to find the offset:
if isempty(find(round(avgAPafterPeak)==0, 1))
    avgAPoffsetInMilliSeconds = xSubset(end);
    avgAPoffsetInDataPoints = length(xSubset); 
end


%% PLOT - AP width ==============================================

% Plot all APs and avg AP (not filtered, baseline subtracted)
% Ddy based ONset is marked with a blue arrow >
% Avg based PEAK is marked with a red arrow ^
% Avg based VALLEY is marked with a red arrow v
% Ddy based OFFset is marked with a blue arrow <
% Avg based OFFset is marked with a red arrow <
figure('name', strcat(fileName, " ", analysisDate, ' - firing_vs_light_test - AP width'));     
hold on;
    plot(xSubset, ySubsetAll,'Color', [0.75, 0.75, 0.75, 0.5], 'LineWidth', 0.2);
    plot(xSubset, avgAP,'Color','black','LineWidth',1.5); 
    plot(xSubset(avgAPpeakInDataPoints), max(avgAP),'^','color', 'r');
    plot(xSubset(avgAPvalleyInDataPoints), min(avgAP),'v','color', 'r');
    plot(ddyBasedOffsetInMilliSeconds, avgAP(ddyBasedOffsetInDataPoints), '<', 'color', 'b');
    plot(ddyBasedOnsetInMilliSeconds, avgAP(ddyBasedOnsetInDataPoints), '>', 'color', 'b');
    plot(avgAPoffsetInMilliSeconds, avgAP(avgAPoffsetInDataPoints), '<', 'color', 'r');
    xlabel('Time (ms)');
    ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
    title([fileName ' AP width'],'Interpreter','none');
hold off;
set(gcf,'Position',[1500 50 400 400]);

% Display y axis inverted to match how extracelullar spikes are most often displayed in the literature
set(gca, 'YDir','reverse');

% adding scale bar
xmaxHere = 1000*(preAPinSeconds+postAPinSeconds);
line([xmaxHere-2 xmaxHere],[min(avgAP) min(avgAP)],'Color','k')
line([xmaxHere xmaxHere],[min(avgAP) min(avgAP)+10],'Color','k')
text(xmaxHere-2, min(avgAP)+5, "2 ms")
text(xmaxHere-2, min(avgAP)+10, strcat("10 ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))


%% PLOT - AP width with derivatives =============================

% Plot the first (blue) and second (red) derivative of the avg
% Use this dor troubleshooting and adjusting ddyPeakThreshold and ddyValleyThreshold
figure('name', strcat(fileName, " ", analysisDate, ' - firing_vs_light_test - AP width ddy'));
hold on;
    plot(xSubset, avgAP,'Color','black','LineWidth',1);
    plot(xForDy, dy,'Color', 'b', 'LineWidth', 1);
    plot(xForDdy, ddy,'Color', 'r', 'LineWidth', 1);
    plot(ddyBasedOnsetInMilliSeconds, ddy(ddyBasedOnsetInDataPoints+1), '>', 'color', 'b');
    plot(ddyBasedOffsetInMilliSeconds, ddy(ddyBasedOffsetInDataPoints), '<', 'color', 'b');
    line([0 xSubset(end)],[ddyPeakThreshold ddyPeakThreshold], 'LineStyle', '-.');
    line([0 xSubset(end)],[-ddyValleyThreshold -ddyValleyThreshold], 'LineStyle', '--');
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title([fileName ' AP width ddy'],'Interpreter','none');
    legend('avg AP', 'avg AP dy', 'avg AP ddy', 'ddy based ONset', 'ddy based OFFset', 'ddyPeakThreshold', 'ddyValleyThreshold', 'Location', 'northeast');
hold off;
set(gcf,'Position',[1500 350 400 400]);


end