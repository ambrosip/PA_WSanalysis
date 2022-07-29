%{ 
DOCUMENTATION
Created: 2022 07 21
Author: PA

This function is used to analyze light-evoked changes in firing rate in
sCRACM experiments using the polygon.
This function was derived from firing_vs_light_ch.

INPUTS explained:
    - gridSize: number of columns in the square grid

    - sweepsPerSquare: number of repeats for each illuminated
    location/square.

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

    - lightChannel: channel where the info about the light stim is stored.
    Usually 2 or 3.

    - singleLightPulse: boolean variable (0 or 1) determining whether this
    sweep has a single light pulse (1) or a train of light pulses (0).
    Added to avoid errors

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
    lightChannel = 2;
    singlelightpulse = 0;

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
    Fig
    XLS

ASSUMPTIONS: 
    - Recording was done in cell attached or loose cell mode in VClamp -
    bandpass filter is appropriate. 
    - For Avg XLS file: light stim parameters are the same in all sweeps.
    - Dopaminergic cells have total AP duration > 2 ms.
    - Irregular cells have ISI CV > 0.2

BEWARE:
    - if you add more variables into "data" for exporting, you need to
    adjust the code for lightEffect.

TO DO:
    - create 'test' version with sweep by sweep plots - DONE
    - complete documentation - ALMOST DONE
%}

% function firing_vs_light_polygon(obj)
obj = m571.s0138;


%%  USER INPUT ==================================================

% Affects data analysis - Organizing data by o-stim grid
gridColumns = 5;
gridRows = 5;
sweepsPerSquare = 3;

% Affects data analysis - Finding APs:
discardedSweeps = [];
discardedSweepsFromEnd = 0;
peaksOrValleys = 'v';   
highpassThreshold = 100;
lowpassThreshold = 1500;    % was 1500
minPeakHeight = 20;         
minPeakDistance = 0.001; 
lightExtensionFactor = 1;
lightChannel = 2;
singleLightPulse = 1; 

% Affects data analysis - AP shape:
preAPinSeconds = 0.005;            
postAPinSeconds = 0.01;           
preAPbaselineDurationSeconds = 0.002;
ddyValleyThreshold = 600;
ddyPeakThreshold = 300;
  
% Affects data display: 
ymax = 200;
ymaxhist = 30;
zoomWindow = 0.25;
ymaxIsiCV = 150;

% Affects data saving:
savefileto = 'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries\2022\2022-07-12 polygon DATs';


%% PREP - get info from file and create arrays ==================

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% getting info from file
mouseNumber = getMouseNumber(obj);
experimentDate = getExperimentDate(obj);
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

% creating matrixes/arrays that will be filled
yFilteredAll = [];
allTimeStamps = [];
baselineTimeStamps = [];
allIsiBaseline = [];
hzBaselineBySweep = [];
hzPreLightBySweep = [];
hzDuringLightBySweep = [];
hzPostLightBySweep = [];
isiMeanBySweep = [];
isiStdBySweep = [];
isiCvBySweep = [];
isIrregularBySweep = [];
data = [];

% creating matrixes/arrays that will be filled for AP shape
xSubset = [];
ySubsetForAPshape = [];
yFilteredSubset = [];
xSubsetAll = [];
ySubsetAll = [];
yFilteredSubsetAll = [];
nAPtotal = 0;

% creating cell arrays that will be filled (can take columns of different sizes)
tsBySweep = {};
isiBySweep = {};
sweepNumberArrayBySweep = {};


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
    [xch2,ych2] = obj.xy(sweepNumber, lightChannel);      
    
    % get light stim parameters
    lightPulseStart = find(diff(ych2>1)>0);
    lightPulseEnd = find(diff(ych2<1)>0);
    lightOnsetTime = lightPulseStart(1)/samplingFrequency;                       % in seconds
    stimDur = (lightPulseEnd(end)-lightPulseStart(end))/samplingFrequency;       % duration of each light pulse in the train (s)

    % if the light stim is a train (singleLightPulse = 0), compute light
    % train information. 
    if singleLightPulse == 0
        stimInterval = (lightPulseStart(2)-lightPulseStart(1))/samplingFrequency;    % interval between each pulse (s)
        stimFreq = 1/stimInterval;                                                   % frequency of the light stim (Hz)
        lightDur = (lightPulseStart(end)-lightPulseStart(1))/samplingFrequency + stimInterval;  % duration of the whole light train stim (s) 
    % if the light stim is a single pulse (singleLightPulse = 1), set the
    % train information to the following values (to avoid errors)
    else
        stimInterval = 0;
        stimFreq = 1;
        lightDur = stimDur;
    end        
    
    % saving data for niceplot
    % y data for each sweep is in a column
    yFilteredAll = [yFilteredAll, yFiltered];    
    
    % create list of sweepNumber with the same size as list of timestamps 
    % to organize raster plot 
    sweepNumberArray = sweepNumber.* ones(length(locs),1);
    
    % storing all AP timestamps
    allTimeStamps = [allTimeStamps; locs];
    
    % getting all of the APs timestamps prior to light stim (aka full
    % baseline, from 0s to lightOnsetTime)
    % to calculate mean baseline Hz and SD later
    locsBaseline = locs;
    indicesToDelete = find(locs >= lightOnsetTime);
    locsBaseline(indicesToDelete) = [];
    baselineTimeStamps = [baselineTimeStamps; locsBaseline];
    
    % getting all the ISIs (inter-spike-intervals) prior to light stim
    % to plot baseline ISI histogram and calculate ISI CV
    isiBaseline = diff(locsBaseline);
    allIsiBaseline = [allIsiBaseline; isiBaseline];
    
%     % storing sweep by sweep data in a structure
%     % fyi access sweep 1 data using dataPerSweep.s0128ts
%     dataPerSweep.(strcat('s',num2str(sweepNumber),'ts')) = locs;
%     dataPerSweep.(strcat('s',num2str(sweepNumber),'hz')) = length(locsBaseline)/lightOnsetTime;
%     dataPerSweep.(strcat('s',num2str(sweepNumber),'isi')) = isiBaseline;
    %----------------------------------------------------------------
    
    % Getting peak/valley amplitudes pre, during and post light
    % to do manual quality control of found peaks later
    if peaksOrValleys == 'peaks'
        pksPreLight = pks;
        pksDuringLight = pks;
        pksPostLight = pks;
        peaksOrValleysAsNum = 1;
    else
        pksPreLight = -pks;
        pksDuringLight = -pks;
        pksPostLight = -pks;
        peaksOrValleysAsNum = -1;
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

    % Data storage    
    % storing sweep by sweep data in a cell array
    % to export later
    % fyi access sweep 1 data using cell2mat(sweepData(1))
    tsBySweep = [tsBySweep, locs];
    isiBySweep = [isiBySweep, isiBaseline];
    sweepNumberArrayBySweep = [sweepNumberArrayBySweep, sweepNumberArray];
    
    % storing sweep by sweep data in an array for easy mean & std calculations later
    hzBaselineBySweep = [hzBaselineBySweep, length(locsBaseline)/lightOnsetTime];
    hzPreLightBySweep = [hzPreLightBySweep, length(locsPreLight)/lightDur];
    hzDuringLightBySweep = [hzDuringLightBySweep, length(locsDuringLight)/lightDur];
    hzPostLightBySweep = [hzPostLightBySweep, length(locsPostLight)/lightDur];
    isiCvBySweep = [isiCvBySweep, std(isiBaseline)/mean(isiBaseline)];
    isiMeanBySweep = [isiMeanBySweep, mean(isiBaseline)];
    isiStdBySweep = [isiStdBySweep, std(isiBaseline)];
    %----------------------------------------------------------------
    
    % checking if cell is irregular (ISI CV > 0.2)
    % ASSUMPTION ALERT, MIGHT NEED UPDATING
    if std(isiBaseline)/mean(isiBaseline) > 0.2
        isIrregular = 1;
    else 
        isIrregular = 0;
    end
    
    isIrregularBySweep = [isIrregularBySweep, isIrregular];
    %----------------------------------------------------------------
    
    % Collecting AP shape data   
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
    
    % Data that will be exported       
    % storing a subset of sweep by sweep data that will be exported
    % this dataset is the one I'm used to manipulating with a few changes
    data = [data; ...
        mouseNumber, ...
        experimentDate, ...
        sweepNumber, ...
        discardedSweepsFromEnd, ...
        peaksOrValleysAsNum, ...
        highpassThreshold, ...
        lowpassThreshold, ...
        minPeakHeight, ...        
        minPeakDistance, ...    
        lightExtensionFactor, ...
        stimDur, ... 
        stimFreq, ...
        lightDur, ...
        length(locsBaseline)/lightOnsetTime, ...
        mean(isiBaseline), ...
        std(isiBaseline), ...
        std(isiBaseline)/mean(isiBaseline), ...
        isIrregular, ...
        length(locsPreLight), ...
        length(locsDuringLight), ...
        length(locsPostLight), ...
        length(locsPreLight)/lightDur, ...
        length(locsDuringLight)/lightDur, ...
        length(locsPostLight)/lightDur];   
    %----------------------------------------------------------------
    
end


%% CELL ANALYSIS - firing (all sweeps, irrespective of square) =======================================

% Mean and Std for pre-light baseline firing rate
hzPreLightMean = mean(hzPreLightBySweep);
hzPreLightStd = std(hzPreLightBySweep);

% counting APs accross ALL SWEEPS
edges = [0:30];
[N, edges] = histcounts(allTimeStamps,edges);
firingHz = N/length(allSweeps);

% is firing modulated by light?
% if cell is inhibited, lightEffect = -1
% if cell is excited, lightEffect = 1
% if cell is indifferent, lightEffect = 0
% data(sweepNumber, 23) is duringLightHz
% data(sweepNumber, 22) is preLightHz
lightEffect = [];
sdFromPreLightHz = [];
for sweepNumber = [1:length(allSweeps)]
    sdFromPreLightHz = [sdFromPreLightHz; (data(sweepNumber, 23) - data(sweepNumber, 22)) / hzPreLightStd];
    if data(sweepNumber, 23) < hzPreLightMean - 2*hzPreLightStd
        lightEffect = [lightEffect; -1];
    elseif data(sweepNumber, 23) > hzPreLightMean + 2*hzPreLightStd
        lightEffect = [lightEffect; +1];
    else 
        lightEffect = [lightEffect; 0];
    end
end

% add lightEffect and sdFromPreLightHz as the last columns of the sweep by sweep data
data = [data, lightEffect, sdFromPreLightHz];




%% CELL ANALYSIS - AP shape (all sweeps, irrespective of square)=====================================

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
%----------------------------------------------------------------

% Calculate AP width and duration based on multiple criteria
halfWidth = avgAPpeakInMilliSeconds - avgAPvalleyInMilliSeconds;
biphasicDuration = avgAPpeakInMilliSeconds - ddyBasedOnsetInMilliSeconds;
totalDurationDdyBased = ddyBasedOffsetInMilliSeconds - ddyBasedOnsetInMilliSeconds;
totalDurationAvgBased = avgAPoffsetInMilliSeconds - ddyBasedOnsetInMilliSeconds;
%----------------------------------------------------------------

% Check if AP duration is consistent with DA cell
% DA cells have total duration > 2 ms
% To be conservative, I'm taking the average between my two duration
% metrics
if (totalDurationDdyBased + totalDurationAvgBased)/2 > 2
    isDA = 1;
else
    isDA = 0;
end
%----------------------------------------------------------------

% store AP shape data 
% cannot add to sweep by sweep data cuz it's just 1 row
dataAPshape = [mouseNumber, ...
    experimentDate, ...
    firstSweepNumber, ...
    lastSweepNumber, ...
    discardedSweepsFromEnd, ...
    peaksOrValleysAsNum, ...
    highpassThreshold, ...
    lowpassThreshold, ...
    minPeakHeight, ...        
    minPeakDistance, ... 
    preAPinSeconds, ...
    postAPinSeconds, ...
    preAPbaselineDurationSeconds, ...
    ddyValleyThreshold, ...
    ddyPeakThreshold, ...
    nAPtotal, ...
    ddyBasedOnsetInMilliSeconds, ...
    avgAPvalleyInMilliSeconds, ...
    avgAPpeakInMilliSeconds, ...
    ddyBasedOffsetInMilliSeconds, ...
    avgAPoffsetInMilliSeconds, ...
    halfWidth, ...
    biphasicDuration, ...
    totalDurationDdyBased, ...
    totalDurationAvgBased, ...
    isDA];


%% CELL ANALYSIS - POLYGON SPECIFIC CODE - avg for each square ==========================

% Assigning sweeps to each o-stim area (square)
% Each row in sweepIDperSquare correspond to a square
% Each column in sweepIDperSquare corresponds to a repeat sweep in that
% square.
totalSquares = gridColumns * gridRows;
sweepOrderPerSquare = []; % from 1 to totalSquares
sweepNumberPerSquare = []; % from 1st sweep to last sweep number
dataSquare = [];

% each row is a square
for row=[1:totalSquares]    
    
    % each column is a sweep repeat for a square
    for column=[1:sweepsPerSquare]    
        
        % assign sweepID 1 to total # of sweeps
        sweepOrderPerSquare(row,column) = row + totalSquares * (column - 1);   
    end
    
    % create a copy of the data per sweep
    dataSubset = data;
    
    % remove all data from sweeps not included in this square
    dataSubset(setdiff(1:end,sweepOrderPerSquare(row,:)),:) = [];
    
    % store average data per square
    dataSquare = [dataSquare; mean(dataSubset)];
end

% adjust sweepIDs so that sweep#1 = (1 + actual number for sweep#1)
sweepNumberPerSquare = sweepOrderPerSquare + allSweeps(1) - 1;  

% concatenate dataSquare with sweepNumberPerSquare to store the number of
% indivudual sweeps that were averaged per square
dataSquare = [dataSquare, sweepNumberPerSquare];






isIrregularCell = median(isIrregularBySweep);
lightEffectCell = median(lightEffect);
sdFromPreLightHzCell = median(sdFromPreLightHz);

dataCell = [mouseNumber, ...
    experimentDate, ...
    firstSweepNumber, ...
    lastSweepNumber, ...
    length(allSweeps), ...
    discardedSweepsFromEnd, ...
    peaksOrValleysAsNum, ...
    highpassThreshold, ...
    lowpassThreshold, ...
    minPeakHeight, ...        
    minPeakDistance, ...    
    lightExtensionFactor, ...
    preAPinSeconds, ...
    postAPinSeconds, ...
    preAPbaselineDurationSeconds, ...
    ddyValleyThreshold, ...
    ddyPeakThreshold, ...    
    stimDur, ... 
    stimFreq, ...
    lightDur, ...    
    nAPtotal, ...
    halfWidth, ...
    biphasicDuration, ...
    totalDurationDdyBased, ...
    totalDurationAvgBased, ...    
    mean(hzBaselineBySweep), ...
    std(hzBaselineBySweep), ...
    std(hzBaselineBySweep)/mean(hzBaselineBySweep), ...
    mean(allIsiBaseline), ...
    std(allIsiBaseline), ...
    std(allIsiBaseline)/mean(allIsiBaseline), ...
    mean(hzPreLightBySweep), ...
    std(hzPreLightBySweep), ...
    mean(hzDuringLightBySweep), ...
    std(hzDuringLightBySweep), ...
    mean(hzPostLightBySweep), ...
    std(hzPostLightBySweep), ...
    (mean(hzDuringLightBySweep) - mean(hzPreLightBySweep))/std(hzPreLightBySweep), ...
    mean(hzDuringLightBySweep)/mean(hzPreLightBySweep), ...
    isDA, ...
    isIrregularCell, ...
    lightEffectCell, ...
    sdFromPreLightHzCell];


% % % %% PLOT - ISI ===================================================
% % % 
% % % % ISI CV pre-light across all sweeps
% % % figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light - baseline ISI counts'));
% % % 
% % % % edges: from 0-1s in 1ms steps
% % % edges = [0:0.001:1]; 
% % % 
% % % histogram(allIsiBaseline, edges);
% % % title([strcat(fileName, ' baseline ISI counts')],'Interpreter','none');
% % % xlabel('ISI (s)');
% % % ylabel('Counts per bin');
% % % axis([0 1 0 ymaxIsiCV])
% % % xticks([0 1]);
% % % set(gcf,'Position',[50 550 400 400]);
% % % yticks([0 ymaxIsiCV]);
% % % 
% % % 
% % % %% PLOT - ISI normalized ========================================
% % % 
% % % % ISI CV pre-light across all sweeps
% % % figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light - baseline ISI prob'));
% % % 
% % % % edges: from 0-1s in 1ms steps
% % % edges = [0:0.001:1]; 
% % % 
% % % histogram(allIsiBaseline, edges, 'Normalization', 'probability', 'FaceColor', [0 0 0],'FaceAlpha', 1, 'EdgeColor', 'none');  
% % % title([strcat(fileName, ' baseline ISI prob')],'Interpreter','none');
% % % xlabel('ISI (s)');
% % % ylabel('Counts ber Bin / Total Counts');
% % % axis([0 0.5 0 0.1])
% % % xticks([0 0.5]);
% % % yticks([0 0.1]);
% % % set(gcf,'Position',[50 50 400 400]);
% % % 
% % % 
% % % % No longer need this since the sweeps are shorter (8s instead of 30s) when
% % % % using a polygon grid
% % % % %% PLOT - Raster plot and histogram =============================
% % % % % Complete - mean +- 2SD Hz is from long baseline (all data prior to light stim)
% % % % 
% % % % figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light - raster and hist'));
% % % % subplot(2,1,1)
% % % % hold on;
% % % % 
% % % % % plotting raster plot 
% % % % for sweep = 1:length(allSweeps)
% % % %     plot(cell2mat(tsBySweep(sweep)), cell2mat(sweepNumberArrayBySweep(sweep)), '|', 'Color', 'k')
% % % % end
% % % % 
% % % % % beautifying raster plot
% % % % title([strcat(fileName, ' raster and hist')],'Interpreter','none');
% % % % axis([0 30 firstSweepNumber-1 lastSweepNumber+1])
% % % % ylabel(strcat('Sweeps (', num2str(length(allSweeps)), ')'));
% % % % yticks([]);
% % % % xticks([]);
% % % % 
% % % % % adding light stim
% % % % rectangle('Position', [lightOnsetTime firstSweepNumber-1 lightDur lastSweepNumber+1], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');
% % % % 
% % % % % stop plotting things on this subplot
% % % % hold off;
% % % % 
% % % % % flip the y-axis so that the first sweep is at the top and the last
% % % % % sweep is at the bottom
% % % % set(gca, 'YDir','reverse');
% % % % 
% % % % % plotting Hz histogram 
% % % % subplot(2,1,2)
% % % % hold on;
% % % % 
% % % % histogram('BinEdges', 0:30, 'BinCounts', firingHz, 'DisplayStyle', 'stairs', 'EdgeColor', 'k'); % 'EdgeColor','none',
% % % % 
% % % % % plot light stim as rectangle
% % % % rectangle('Position', [lightOnsetTime 0 lightDur ymaxhist], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');
% % % % 
% % % % % plot Hz mean as horizontal line
% % % % yline(hzBaselineMean, '--');
% % % % 
% % % % % plot +- 2 SD as rectangle around mean
% % % % % [x y width height]
% % % % rectangle('Position', [0 hzBaselineMean-(2*hzBaselineStd) 30 4*hzBaselineStd], 'FaceColor', [0 0 0 0.1], 'EdgeColor', 'none');
% % % % 
% % % % xlabel('Time (s)');
% % % % ylabel('Firing rate (Hz)');
% % % % axis([0 30 0 ymaxhist])
% % % % set(gcf,'Position',[500 550 500 400]);
% % % % % xticks([0 30]);
% % % % yticks([0 ymaxhist]);
% % % % hold off;
% % % 
% % % 
% % % %% PLOT - Raster plot and histogram =============================
% % % % Zoomed in - mean +- 2SD Hz is from short pre-light baseline
% % % 
% % % figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light - raster and hist zoom'));
% % % subplot(2,1,1)
% % % hold on;
% % % 
% % % % plotting raster plot 
% % % for sweep = 1:length(allSweeps)
% % %     plot(cell2mat(tsBySweep(sweep)), cell2mat(sweepNumberArrayBySweep(sweep)), '|', 'Color', 'k')
% % % end
% % % 
% % % % zooming in and beautifying raster plot
% % % title([strcat(fileName, ' raster and hist zoom')],'Interpreter','none');
% % % axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur firstSweepNumber-1 lastSweepNumber+1])
% % % ylabel(strcat('Sweeps (', num2str(length(allSweeps)), ')'));
% % % yticks([]);
% % % xticks([]);
% % % 
% % % % adding light stim
% % % rectangle('Position', [lightOnsetTime firstSweepNumber-1 lightDur lastSweepNumber+1], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');
% % % 
% % % % stop plotting things on this subplot
% % % hold off;
% % % 
% % % % flip the y-axis so that the first sweep is at the top and the last
% % % % sweep is at the bottom
% % % set(gca, 'YDir','reverse');
% % % 
% % % % plotting Hz histogram 
% % % subplot(2,1,2)
% % % hold on;
% % % 
% % % histogram('BinEdges', 0:30, 'BinCounts', firingHz, 'DisplayStyle', 'stairs', 'EdgeColor', 'k'); % 'EdgeColor','none',
% % % 
% % % % plot light stim as rectangle
% % % rectangle('Position', [lightOnsetTime 0 lightDur ymaxhist], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');
% % % 
% % % % plot Hz mean as horizontal line
% % % yline(hzPreLightMean, '--');
% % % 
% % % % plot +- 2 SD as rectangle around mean
% % % % [x y width height]
% % % rectangle('Position', [0 hzPreLightMean-(2*hzPreLightStd) 30 4*hzPreLightStd], 'FaceColor', [0 0 0 0.1], 'EdgeColor', 'none');
% % % 
% % % xlabel('Time (s)');
% % % ylabel('Firing rate (Hz)');
% % % axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur 0 ymaxhist])
% % % set(gcf,'Position',[500 50 500 400]);
% % % %     xticks([0 30]);
% % % yticks([0 ymaxhist]);
% % % hold off;
% % % % movegui('east')
% % % 
% % % 
% % % %% PLOT - AP width ==============================================
% % % 
% % % % Plot all APs and avg AP (not filtered, baseline subtracted)
% % % % Ddy based ONset is marked with a blue arrow >
% % % % Avg based PEAK is marked with a red arrow ^
% % % % Avg based VALLEY is marked with a red arrow v
% % % % Ddy based OFFset is marked with a blue arrow <
% % % % Avg based OFFset is marked with a red arrow <
% % % figure('name', strcat(fileName, " ", analysisDate, ' - firing_vs_light - AP width'));     
% % % hold on;
% % %     plot(xSubset, ySubsetAll,'Color', [0.75, 0.75, 0.75, 0.5], 'LineWidth', 0.2);
% % %     plot(xSubset, avgAP,'Color','black','LineWidth',1.5); 
% % %     plot(xSubset(avgAPpeakInDataPoints), max(avgAP),'^','color', 'r');
% % %     plot(xSubset(avgAPvalleyInDataPoints), min(avgAP),'v','color', 'r');
% % %     plot(ddyBasedOffsetInMilliSeconds, avgAP(ddyBasedOffsetInDataPoints), '<', 'color', 'b');
% % %     plot(ddyBasedOnsetInMilliSeconds, avgAP(ddyBasedOnsetInDataPoints), '>', 'color', 'b');
% % %     plot(avgAPoffsetInMilliSeconds, avgAP(avgAPoffsetInDataPoints), '<', 'color', 'r');
% % %     xlabel('Time (ms)');
% % %     ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
% % %     title([fileName ' AP width'],'Interpreter','none');
% % % hold off;
% % % set(gcf,'Position',[1000 50 400 400]);
% % % 
% % % % Display y axis inverted to match how extracelullar spikes are most often displayed in the literature
% % % set(gca, 'YDir','reverse');
% % % 
% % % % adding scale bar
% % % xmaxHere = 1000*(preAPinSeconds+postAPinSeconds);
% % % line([xmaxHere-2 xmaxHere],[min(avgAP) min(avgAP)],'Color','k')
% % % line([xmaxHere xmaxHere],[min(avgAP) min(avgAP)+10],'Color','k')
% % % text(xmaxHere-2, min(avgAP)+5, "2 ms")
% % % text(xmaxHere-2, min(avgAP)+10, strcat("10 ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
% % % 
% % % 
% % % %% PLOT - AP width with derivatives =============================
% % % 
% % % % Plot the first (blue) and second (red) derivative of the avg
% % % % Use this dor troubleshooting and adjusting ddyPeakThreshold and ddyValleyThreshold
% % % figure('name', strcat(fileName, " ", analysisDate, ' - firing_vs_light - AP width ddy'));
% % % hold on;
% % %     plot(xSubset, avgAP,'Color','black','LineWidth',1);
% % %     plot(xForDy, dy,'Color', 'b', 'LineWidth', 1);
% % %     plot(xForDdy, ddy,'Color', 'r', 'LineWidth', 1);
% % %     plot(ddyBasedOnsetInMilliSeconds, ddy(ddyBasedOnsetInDataPoints+1), '>', 'color', 'b');
% % %     plot(ddyBasedOffsetInMilliSeconds, ddy(ddyBasedOffsetInDataPoints), '<', 'color', 'b');
% % %     line([0 xSubset(end)],[ddyPeakThreshold ddyPeakThreshold], 'LineStyle', '-.');
% % %     line([0 xSubset(end)],[-ddyValleyThreshold -ddyValleyThreshold], 'LineStyle', '--');
% % %     xlabel('Time (ms)');
% % %     ylabel('Amplitude');
% % %     title([fileName ' AP width ddy'],'Interpreter','none');
% % %     legend('avg AP', 'avg AP dy', 'avg AP ddy', 'ddy based ONset', 'ddy based OFFset', 'ddyPeakThreshold', 'ddyValleyThreshold', 'Location', 'northeast');
% % % hold off;
% % % set(gcf,'Position',[1000 550 400 400]);
% % % 
% % % 
% % % %% PLOT - AP width with derivatives normalized ==================
% % % 
% % % % % Plot the first (blue) and second (red) derivative of the avg - Normalized
% % % % figure('name', strcat(fileName,' y dy ddy Normalized'));
% % % % hold on;
% % % %     plot(xSubset, avgAP/max(avgAP),'Color','black','LineWidth',1);
% % % %     plot(xForDy, dy/max(dy),'Color', 'b', 'LineWidth', 1);
% % % %     plot(xForDdy, ddy/max(ddy),'Color', 'r', 'LineWidth', 1);
% % % %     plot(ddyAfterLastPeakOrValleyInMilliSeconds, ddy(ddyAfterLastPeakOrValleyInDataPoints)/max(ddy), 'o', 'color', 'k');
% % % %     plot(ddyBasedOffsetInMilliSeconds, ddy(ddyBasedOffsetInDataPoints)/max(ddy), '<', 'color', 'b');
% % % %     plot(ddyBasedOnsetInMilliSeconds, ddy(ddyBasedOnsetInDataPoints)/max(ddy), '>', 'color', 'b');
% % % %     xlabel('Time (ms)');
% % % %     ylabel('Normalized Amplitude (au)');
% % % %     title([fileName ' Normalized y dy ddy'],'Interpreter','none');
% % % %     legend('avg AP', 'avg AP dy', 'avg AP ddy', 'Location', 'northwest');
% % % % hold off;
% % % % movegui('northeast');
% % % 
% % % 
% % % %% PLOT - niceplot of first sweep ===============================
% % % 
% % % % plotting niceplot of first sweep    
% % % figure('name', strcat(fileName, " ", analysisDate, ' - firing_vs_light - niceplot 1st sweep'));
% % % plot(x,yFilteredAll(:,1),'k','LineWidth',1);
% % % axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur -ymax ymax]);
% % % title([fileName ' - firing_vs_light - niceplot 1st sweep'],'Interpreter','none');
% % % set(gca,'Visible','off');
% % % set(gcf,'Position',[1400 50 500 400]);
% % % 
% % % % if the light stim is a train (singleLightPulse = 0), add light
% % % % stim cartoon as a train.
% % % % adding light stim as train - faithful to light stim parameters
% % % % note that this code will use light stim parameters from the last sweep!
% % % % if light stim is not the same accross all sweeps, this will be
% % % % misleading!
% % % if singleLightPulse == 0
% % %     for nStim=1:length(lightPulseStart)
% % %         line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[ymax,ymax],'Color',[0 0.4470 0.7410],'LineWidth',10)
% % %     end
% % % end
% % % 
% % % % adding light stim as rectangle
% % % % note that this code will use light stim parameters from the last sweep!
% % % % if light stim is not the same accross all sweeps, this will be
% % % % misleading!
% % % line([lightOnsetTime,lightOnsetTime+lightDur],[ymax-10,ymax-10],'Color',[0 0.4470 0.7410],'LineWidth',10)
% % % 
% % % % adding scale bar
% % % xmin = lightOnsetTime-lightDur;
% % % xmax = lightOnsetTime+2*lightDur;
% % % ymin = -ymax;
% % % line([xmax-1 xmax],[ymax ymax],'Color','k')
% % % line([xmax xmax],[ymax ymax-((ymax-ymin)/10)],'Color','k')
% % % text(xmax-1, ymax-((ymax-ymin)/20), "1 s")
% % % text(xmax-1, ymax-((ymax-ymin)/10), strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
% % % hold off;
% % % 
% % % 
% % % %% PLOT - niceplot of all sweeps ================================
% % % 
% % % % plotting niceplot     
% % % figure('name', strcat(fileName, " ", analysisDate, ' - firing_vs_light - niceplot all'));
% % % plot(x,yFilteredAll,'Color',[0, 0, 0, 0.25]);
% % % axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur -ymax ymax]);
% % % title([fileName ' - firing_vs_light - niceplot all'],'Interpreter','none');
% % % set(gca,'Visible','off');
% % % set(gcf,'Position',[1400 550 500 400]);
% % % 
% % % % adding light stim   
% % % % if the light stim is a train (singleLightPulse = 0), add light
% % % % stim cartoon as a train.
% % % % note that this code will use light stim parameters from the last sweep!
% % % % if light stim is not the same accross all sweeps, this will be
% % % % misleading!
% % % if singleLightPulse == 0
% % %     for nStim=1:length(lightPulseStart)
% % %         line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[ymax,ymax],'Color',[0 0.4470 0.7410],'LineWidth',10)
% % %     end
% % % end
% % % 
% % % % adding light stim
% % % % note that this code will use light stim parameters from the last sweep!
% % % % if light stim is not the same accross all sweeps, this will be
% % % % misleading!
% % % line([lightOnsetTime,lightOnsetTime+lightDur],[ymax-10,ymax-10],'Color',[0 0.4470 0.7410],'LineWidth',10)
% % % 
% % % % adding scale bar
% % % xmin = lightOnsetTime-lightDur;
% % % xmax = lightOnsetTime+2*lightDur;
% % % ymin = -ymax;
% % % line([xmax-1 xmax],[ymax ymax],'Color','k')
% % % line([xmax xmax],[ymax ymax-((ymax-ymin)/10)],'Color','k')
% % % text(xmax-1, ymax-((ymax-ymin)/20), "1 s")
% % % text(xmax-1, ymax-((ymax-ymin)/10), strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
% % % 
% % % 
% % % %% EXPORTING XLS files ==========================================
% % % 
% % % % stores key sweep by sweep data
% % % filename = strcat(fileName, '_', analysisDate, " - firing_vs_light - sweep_by_sweep");
% % % fulldirectory = strcat(savefileto,'\',filename,'.xls');        
% % % dataInCellFormat = {};
% % % dataInCellFormat = num2cell(data);
% % % labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
% % %     {'mouse', ...
% % %     'date', ...
% % %     'sweep', ...
% % %     'discardedSweeps', ...
% % %     'peaks(1)OrValleys(-1)', ...
% % %     'highpassThreshold', ...
% % %     'lowpassThreshold', ...
% % %     'minPeakHeight', ...        
% % %     'minPeakDistance', ...    
% % %     'lightExtensionFactor', ...
% % %     'lightPulseDur(s)', ... 
% % %     'lightStimFreq(Hz)', ...
% % %     'lightDur(s)', ...
% % %     'baselineHz', ...
% % %     'baselineIsiMean(s)', ...
% % %     'baselineIsiSD(s)', ...
% % %     'baselineIsiCV', ...
% % %     'irregular', ...
% % %     'preLightAPs', ...
% % %     'duringLightAPs', ...
% % %     'postLightAPs', ...
% % %     'preLightHz', ...
% % %     'duringLightHz', ...
% % %     'postLightHz', ...
% % %     'lightEffect', ...
% % %     'sdFromPreLightHz'});
% % % writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
% % % disp('I saved the sweep_by_sweep xls file')
% % % 
% % % % stores all timestamps from all sweeps in a single column
% % % filename = strcat(fileName, '_', analysisDate, " - firing_vs_light - all_AP_timestamps");
% % % fulldirectory = strcat(savefileto,'\',filename,'.xls');        
% % % writematrix(allTimeStamps, fulldirectory, 'WriteMode', 'overwritesheet');
% % % disp('I saved the all_AP_timestamps xls file')
% % % 
% % % % stores AP shape data in a single row
% % % filename = strcat(fileName, '_', analysisDate, " - firing_vs_light - AP_shape");
% % % fulldirectory = strcat(savefileto,'\',filename,'.xls');        
% % % dataInCellFormat = {};
% % % dataInCellFormat = num2cell(dataAPshape);
% % % labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
% % %     {'mouse', ...
% % %     'date', ...
% % %     'firstSweep', ...
% % %     'lastSweep', ...
% % %     'discardedSweepsFromEnd', ...
% % %     'peaks(1)OrValleys(-1)', ...
% % %     'highpassThreshold', ...
% % %     'lowpassThreshold', ...
% % %     'minPeakHeight', ...        
% % %     'minPeakDistance', ...
% % %     'preAP(s)', ...
% % %     'postAP(s)', ...
% % %     'baselinePreAPdur(s)', ...
% % %     'ddyValleyThreshold', ...
% % %     'ddyPeakThreshold', ...
% % %     'nAPtotal', ...
% % %     'ddyBasedOnset(ms)', ...
% % %     'APvalley(ms)', ...
% % %     'APpeak(ms)', ...
% % %     'ddyBasedOffset(ms)', ...
% % %     'avgBasedOffset(ms)', ...
% % %     'halfWidth(ms)', ...
% % %     'biphasicDuration(ms)', ...
% % %     'totalDurationDdyBased(ms)', ...
% % %     'totalDurationAvgBased(ms)', ...
% % %     'isDA'});
% % % writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
% % % disp('I saved the AP_shape xls file')
% % % 
% % % % stores avg data in a single row
% % % filename = strcat(fileName, '_', analysisDate, " - firing_vs_light - cell_avgs");
% % % fulldirectory = strcat(savefileto,'\',filename,'.xls');        
% % % dataInCellFormat = {};
% % % dataInCellFormat = num2cell(dataCell);
% % % labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
% % %     {'mouse', ...
% % %     'date', ...
% % %     'firstSweep', ...
% % %     'lastSweep', ...
% % %     'nSweeps', ...
% % %     'discardedSweepsFromEnd', ...
% % %     'peaks(1)OrValleys(-1)', ...
% % %     'highpassThreshold', ...
% % %     'lowpassThreshold', ...
% % %     'minPeakHeight', ...        
% % %     'minPeakDistance', ...
% % %     'lightExtensionFactor', ...
% % %     'preAP(s)', ...
% % %     'postAP(s)', ...
% % %     'baselinePreAPdur(s)', ...
% % %     'ddyValleyThreshold', ...
% % %     'ddyPeakThreshold', ...
% % %     'lightPulseDur(s)', ... 
% % %     'lightStimFreq(Hz)', ...
% % %     'lightDur(s)', ...   
% % %     'nAPtotalForShape', ...
% % %     'halfWidth(ms)', ...
% % %     'biphasicDuration(ms)', ...
% % %     'totalDurationDdyBased(ms)', ...
% % %     'totalDurationAvgBased(ms)', ...
% % %     'baselineHzMean', ...
% % %     'baselineHzSD', ...
% % %     'baselineHzCV', ...
% % %     'baselineISImean', ...
% % %     'baselineISIsd', ...
% % %     'baselineISIcv', ...
% % %     'preLightHzMean', ...
% % %     'preLightHzSD', ...
% % %     'duringLightHzMean', ...
% % %     'duringLightHzSD', ...
% % %     'postLightHzMean', ...
% % %     'postLightHzSD', ...
% % %     'lightEffect(SDfromMean)', ...
% % %     'lightEffect(duringHz/preHz)', ...
% % %     'isDA(0,1)', ...
% % %     'isIrregular(0,1)', ...
% % %     'lightEffect(-1,0,1)', ...
% % %     'medianSdFromPreLightHz'});
% % % writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
% % % disp('I saved the cell_avgs xls file')
% % % 
% % % % print stuff
% % % [isDA, isIrregularCell, lightEffectCell]