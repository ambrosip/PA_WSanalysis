%{ 
DOCUMENTATION
Created: 2021 03 18
Author: PA

This function is used to test the input parameters to be used in the 
function "firing_vs_light".

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

    - lightExtensionFactor: to account for lingering effects after the end
    of the light pulse, extend the time interval considered under "light
    effect". To NOT extend the light pulse (my default), set this variable
    to 1.

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
    lightExtensionFactor = 1;

    preAPinSeconds = 0.005;            
    postAPinSeconds = 0.01;           
    preAPbaselineDurationSeconds = 0.002;
    ddyValleyThreshold = 30;
    ddyPeakThreshold = 10;

    ymax = 75;
    ymaxhist = 15;
    ZoomWindow = 0.25;
    ymaxIsiCV = 150;
    savefileto = 'D:\CORONAVIRUS DATA\From MATLAB';
    
OUTPUTS:
    Figs

ASSUMPTIONS: 
    - Recording was done in cell attached or loose cell mode in VClamp -
    bandpass filter is appropriate. 

TO DO:
    - complete documentation
%}

function firing_vs_light_test(obj)


%%  USER INPUT ==================================================

% Affects data analysis - Finding APs:
discardedSweeps = [];
discardedSweepsFromEnd = 0;
peaksOrValleys = 'v';   
highpassThreshold = 100;
lowpassThreshold = 1500;    
minPeakHeight = 10;         
minPeakDistance = 0.025;    
lightExtensionFactor = 1;

% Affects data analysis - AP shape:
preAPinSeconds = 0.005;            
postAPinSeconds = 0.01;           
preAPbaselineDurationSeconds = 0.002;
ddyValleyThreshold = 60;
ddyPeakThreshold = 35;
  
% Affects data display: 
ymax = 75;
ymaxhist = 15;
zoomWindow = 0.25;
ymaxIsiCV = 150;

% Affects data saving:
savefileto = 'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries\2022\2022-07-12 polygon DATs';


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

    % get light stim data
    [xch2,ych2] = obj.xy(sweepNumber, 2);      
    
    % get light stim parameters
    lightPulseStart = find(diff(ych2>1)>0);
    lightPulseEnd = find(diff(ych2<1)>0);
    lightOnsetTime = lightPulseStart(1)/samplingFrequency;                       % in seconds
    stimDur = (lightPulseEnd(end)-lightPulseStart(end))/samplingFrequency;       % duration of each light pulse in the train (s)
    stimInterval = (lightPulseStart(2)-lightPulseStart(1))/samplingFrequency;    % interval between each pulse (s)
    stimFreq = 1/stimInterval;                                                   % frequency of the light stim (Hz)
    lightDur = (lightPulseStart(end)-lightPulseStart(1))/samplingFrequency + stimInterval;  % duration of the whole light train stim (s)        
    %----------------------------------------------------------------
    
    
    % Getting peak/valley amplitudes pre, during and post light
    % to do manual quality control of found peaks later
    if peaksOrValleys == 'peaks'
        pksPreLight = pks;
        pksDuringLight = pks;
        pksPostLight = pks;
    else
        pksPreLight = -pks;
        pksDuringLight = -pks;
        pksPostLight = -pks;
    end
    %----------------------------------------------------------------
    
    
    % Getting timestamps pre, during, and post light  
        % locsPreLight will be different from locsBaseline: locsPreLight is
        % the baseline immediately before the light stim (from
        % lightOnsetTime-lightDur to lightOnsetTime), while locsBaseline is
        % the full baseline (from 0s to lightOnsetTime).
    % first, store locs (timestamps) in new variables 
    locsPreLight = locs;
    locsDuringLight = locs;
    locsPostLight = locs;
    
    % second, find indices at which the value of locs is smaller than the
    % light onset or larger than the light offset (with optional extension
    % factor to look for lingering light effect).
    indicesToDelete = find(locs<lightOnsetTime | locs>(lightOnsetTime+lightDur*lightExtensionFactor));
    
    % then, delete all peaks found before or after the light pulse
    locsDuringLight(indicesToDelete) = [];
    pksDuringLight(indicesToDelete) = [];
    
    % third, find indices for timestamps beyond the immediate pre-light
    % period and delete unwanted peaks
    indicesToDelete = find(locs<lightOnsetTime-lightDur | locs>=lightOnsetTime);
    locsPreLight(indicesToDelete) = [];
    pksPreLight(indicesToDelete) = [];

    % fourth, find indices for timestamps beyond the immediate post-light
    % period and delete unwanted peaks
    indicesToDelete = find(locs<(lightOnsetTime+lightDur*lightExtensionFactor) | locs>(lightOnsetTime+(lightDur*lightExtensionFactor)+lightDur));
    locsPostLight(indicesToDelete) = [];
    pksPostLight(indicesToDelete) = [];
    %----------------------------------------------------------------
    
    
    % Collecting AP shape data  
    % getting all of the APs timestamps prior to light stim (aka full
    % baseline, from 0s to lightOnsetTime)
    locsBaseline = locs;
    indicesToDelete = find(locs >= lightOnsetTime);
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
    figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_test - ISI histogram'));      
    
        % Ch1 (recorded data) histogram
        subplot(2,1,1)         
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
        subplot(2,1,2)   
        plot(xch2,ych2);
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([0 30 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));        
        set(gcf,'Position',[50 550 280 420]) % top left corner
    %----------------------------------------------------------------        
        
    
    % PLOT firing histogram for each sweep (#APs)
    sweepDuration = obj.header.Acquisition.Duration;
    sweepTime=0;
    nAPperSecBin=[];
    
    % naming figure file    
    figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_test - nAPs histogram'));      
    
        % Ch1 (recorded data) histogram
        subplot(2,1,1)         
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
        subplot(2,1,2)   
        plot(xch2,ych2);
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([0 30 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));        
        set(gcf,'Position',[350 550 280 420])
    %----------------------------------------------------------------  
        
    
    % PLOT found peaks for each sweep    
    hzPreLight = length(locsPreLight)/lightDur;
    hzDuringLight = length(locsDuringLight)/lightDur;
    hzPostLight = length(locsPostLight)/lightDur;   
    
    figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_test - found peaks')); 
        
        % Ch1 subplot (recorded data)
        subplot(2,1,1);
        plot(x,yFiltered)
        hold on;
        plot(locsDuringLight,pksDuringLight,'o','color','red');
        plot(locsPreLight,pksPreLight,'o','color','blue');
        plot(locsPostLight,pksPostLight,'o','color','blue');
        hold off;
        axis([lightOnsetTime-4 lightOnsetTime+lightDur+4 -inf ymax])
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([obj.file ' (' num2str(sweepNumber) ') - peaks'],'Interpreter','none');
        text(lightOnsetTime + lightDur/4,0.9*ymax,num2str(round(hzDuringLight,2)),'color','red');
        text(lightOnsetTime - 3,0.9*ymax,num2str(round(hzPreLight,2)),'color','blue');
        text(lightOnsetTime + lightDur + 2,0.9*ymax,num2str(round(hzPostLight,2)),'color','blue');
        
        % Ch2 subplot (light stim command)
        subplot(2,1,2);
        plot(xch2,ych2)
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([lightOnsetTime-4 lightOnsetTime+lightDur+4 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));
        set(gcf,'Position',[650 550 500 400])
    %----------------------------------------------------------------         
    
    
    % PLOT zoomed in niceplot (start of light stim)
    figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_test - zoom start')); 
    
        % main plot
        plot(x,yFiltered,'k','LineWidth',1)
        xmin = lightOnsetTime-zoomWindow;
        xmax = lightOnsetTime+zoomWindow;
        ymin = -ymax;
        axis([xmin xmax ymin ymax]);
        title([obj.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
        set(gca,'Visible','off');
        set(gcf,'Position',[50 50 500 400])
        
        % adding light stim
        interStimInterval = 1/stimFreq;
        postStimInterval = interStimInterval - stimDur;
        for nStim = 1 : stimFreq * xmax
            startStim = lightOnsetTime + (nStim - 1) * interStimInterval;
            line([startStim, startStim + stimDur],[ymax,ymax],'Color',[0 0.4470 0.7410],'LineWidth',10)
        end
        
        % adding scale bar
        line([xmin,xmin+zoomWindow],[ymin,ymin],'Color','k')
        line([xmin,xmin],[ymin,ymin+((ymax-ymin)/10)],'Color','k')
        text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/30),strcat(num2str(zoomWindow*1000)," ms"))
        text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/12),strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
    %---------------------------------------------------------------- 
     
    
    % PLOT zoomed in niceplot (end of light stim) 
    figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_test - zoom end')); 
    
        % main plot
        plot(x,yFiltered,'k','LineWidth',1)
        xmin=lightOnsetTime+lightDur-zoomWindow;
        xmax=lightOnsetTime+lightDur+zoomWindow;
        ymin=-ymax;
        axis([xmin xmax ymin ymax]);
        title([obj.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
        set(gca,'Visible','off');
        set(gcf,'Position',[600 50 500 400])
        
        % adding light stim
        interStimInterval = 1/stimFreq;
        postStimInterval = interStimInterval-stimDur;
        for nStim=1:stimFreq*zoomWindow
            startStim = lightOnsetTime+lightDur-zoomWindow + (nStim - 1) * interStimInterval;
            line([startStim,startStim+stimDur],[ymax,ymax],'Color',[0 0.4470 0.7410],'LineWidth',10)
        end
        
        % adding scale bar
        line([xmin,xmin+zoomWindow],[ymin,ymin],'Color','k')
        line([xmin,xmin],[ymin,ymin+((ymax-ymin)/10)],'Color','k')
        text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/30),strcat(num2str(zoomWindow*1000)," ms"))
        text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/12),strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
    %----------------------------------------------------------------        
    
    
    % PLOT one figure per sweep that shows whole sweep (raw)    
    figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_test - raw'));   
        
        subplot(2,1,1)
        plot(x,y);
        axis([-inf inf -inf inf]);           
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([' (' num2str(sweepNumber) ') - raw'],'Interpreter','none');
      
        subplot(2,1,2) 
        plot(xch2,ych2);
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([0 30 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));        
        set(gcf,'Position',[1200 550 280 420])
    %----------------------------------------------------------------    

    
    % PLOT one figure per sweep that shows whole sweep (filtered)    
    figure('name', strcat(fileName, ' (', num2str(sweepNumber), ')_', analysisDate, ' - firing_vs_light_test - filtered'));  
        
        subplot(2,1,1)
        plot(x,yFiltered);
        axis([-inf inf -inf ymax]);           
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([' (' num2str(sweepNumber) ') - filtered'],'Interpreter','none');
      
        subplot(2,1,2) 
        plot(xch2,ych2);
        yminhere = min(ych2)-5;
        ymaxhere = max(ych2)+5;
        axis([0 30 yminhere ymaxhere])
        xlabel('Time (s)');
        ylabel(strcat(obj.header.Acquisition.ActiveChannelNames(2), ' (', obj.header.Acquisition.AnalogChannelUnits(2), ')'));        
        set(gcf,'Position',[1200 50 280 420])
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
set(gcf,'Position',[1500 550 400 400]);


end