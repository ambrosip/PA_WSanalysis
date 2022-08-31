%{ 
DOCUMENTATION
Created: 2022 07 21
Author: PA

This function is used to analyze light-evoked changes in firing rate in
sCRACM experiments using the polygon.
This function was derived from firing_vs_light_ch, but went through A LOT
of rewritting.

INPUTS explained:
    - gridSize: number of columns in the square grid

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

    - heatmapMin: sets min value for color-coding the heatmap

    - heatmapMax: sets max value for color-coding the heatmap

    - cellImageFileName: 

    - cellImageDir: folder where cell image file is located

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

    ymax = 75;
    ymaxhist = 15;
    ZoomWindow = 0.25;
    ymaxIsiCV = 150;
    heatmapMin = -2;
    heatmapMax = 0;
    cellImageFileName = 's2c1_dic2.tif';
    cellImageDir = 'E:\Priscilla - BACKUP 20200319\Ephys\2022\20220720 m571 dat nphr';

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
    - An ordered polygon grid design was used - aka not random
    - Width and frequency of polygon stim is the the same as LED o-stim
    (only relevant if you use the polygon channel as the light stim
    channel)

BEWARE:
    - if you add more variables into "data" for exporting, you need to
    adjust the code at multiple places. i think I finally fixed this

TO DO:
    - create 'test' version with sweep by sweep plots - DONE
    - complete documentation - ALMOST DONE
%}

% function firing_vs_light_polygon(obj)
obj = m571.s0138;


%%  USER INPUT ============================================================

% Affects data analysis - Organizing data by o-stim grid
gridColumns = 5;
gridRows = 5;

% Affects data analysis - Finding APs:
discardedSweeps = [157,158];
discardedSweepsFromEnd = 0;
peaksOrValleys = 'v';   
highpassThreshold = 100;
lowpassThreshold = 1500;    
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
  
% Affects data display: 
ymax = 200;
ymaxhist = 15;
zoomWindow = 0.25;
ymaxIsiCV = 150;
heatmapMin = -2;
heatmapMax = 0;
cellImageFileName = 's2c1_dic2.tif';
cellImageDir = 'E:\Priscilla - BACKUP 20200319\Ephys\2022\20220720 m571 dat nphr';

% Affects data saving:
savefileto = 'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries\2022\2022-08-11 polygon';


%% PREP - get monitor info for plot display organization =====================================================

monitorPositions = get(0, 'MonitorPositions' );
maxHeight = monitorPositions(1,4) - 100;
maxWidth = monitorPositions(1,3) - 100;


%% PREP - get info from h5 file and create arrays ============================================================

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% getting info from h5 file
mouseNumber = getMouseNumber(obj);
experimentDate = getExperimentDate(obj);
samplingFrequency = obj.header.Acquisition.SampleRate;
fileName = obj.file;
sweepDurationInSec = obj.header.SweepDuration;
sweepDurationInDataPts = sweepDurationInSec * samplingFrequency;

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

% now checking if the total number of sweeps is a multiple of the total
% number of squares. If not, remove the extra sweeps from the end.
totalSweeps = length(allSweeps);
totalSquares = gridColumns * gridRows;
sweepsPerSquare = totalSweeps / totalSquares;
if mod(totalSweeps, totalSquares) ~= 0
    lastSweepNumber = lastSweepNumber - mod(totalSweeps, totalSquares);
    allSweeps = firstSweepNumber:lastSweepNumber;
    totalSweeps = length(allSweeps);
    sweepsPerSquare = totalSweeps / totalSquares;
end   

% creating matrixes/arrays that will be filled
yFilteredAll = [];
allTimeStamps = [];
baselineTimeStamps = [];
allBaselineIsi = [];
hzBaselineBySweep = [];
hzPreLightBySweep = [];
hzDuringLightBySweep = [];
hzPostLightBySweep = [];
isiMeanBySweep = [];
isiStdBySweep = [];
isiCvBySweep = [];
isIrregularBySweep = [];
data = [];
isDiscardedBySweep = [];

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


%% SWEEP BY SWEEP ANALYSIS ===================================================================================

% get data from all sweeps in file
for sweepNumber = allSweeps
    
    % ignore discarded sweeps
    if ismember(sweepNumber, discardedSweeps)        
        isDiscarded = 1;
        yFiltered = NaN(sweepDurationInDataPts, 1);
        pks = NaN;
        locs = NaN;
            
    % analyze all other sweeps
    else         
        isDiscarded = 0;
        
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
    end
    
    % store info about which sweeps were discarded
    isDiscardedBySweep = [isDiscardedBySweep; isDiscarded];
    %----------------------------------------------------------------------
    
    
    % get light stim data
    [xch2,ych2] = obj.xy(sweepNumber, lightChannel); 

    % get light stim parameters
    % if you are using a TTL pulse (5V) to control the LED, use code #1
    % if you are using an analog output (0-5V), use code #2
    % ALERT let's just use the polygon TTL pulse for now, since the analog
    % output is soooo small, making the detection of the light pulse start
    % and end times really hard with simple methods.
    
    % code #1 - works
    % look for a big change 
    lightPulseStart = find(diff(ych2>1)>0);
    lightPulseEnd = find(diff(ych2<1)>0);

%     % code #2 - does not always work
%     % look for a small change
%     % to avoid artifacts, use both the derivative and the absolute value of
%     % ych2 to find the start and end times of each light pulse
%     % ALERT need to check this code with a train o-stim
%     lightPulseStart = find(diff(ych2)>0.075 & ych2(1:end-1)<0.05);
%     lightPulseEnd = find(diff(ych2)<-0.075 & ych2(1:end-1)>0.05);
    
    % continue to get light stim info
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
    %----------------------------------------------------------------------

    
    % save data for niceplot
    % y data for each sweep is in a column
    yFilteredAll = [yFilteredAll, yFiltered];
    %----------------------------------------------------------------------

    
    % storing all AP timestamps
    allTimeStamps = [allTimeStamps; locs];   

    % getting all of the APs timestamps prior to light stim (aka full
    % baseline, from 0s to lightOnsetTime)
    % to calculate mean baseline Hz and SD later
    locsBaseline = locs;
    indicesToDelete = find(locs >= lightOnsetTime);
    locsBaseline(indicesToDelete) = [];
    baselineTimeStamps = [baselineTimeStamps; locsBaseline];
    %----------------------------------------------------------------------

    
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
     

    % check if there are any NaNs in locs - if so, this is data from a
    % discarded sweep, and we will signal that by assigning NaNs to all
    % stored variables.
    % among the stored variables we have all the ISIs (inter-spike-intervals)
    % prior to light stim.
    % These ISIs will be used to plot the ISI histogram and calculate ISI CV 
    if any(isnan(locs))       
        hzBaseline = NaN;
        hzPreLight = NaN;
        hzDuringLight = NaN;
        hzPostLight = NaN;
        isiBaseline = NaN;
        isiCv = NaN;
        isiMean = NaN;
        isiStd = NaN;        
    else        
        hzBaseline = length(locsBaseline)/lightOnsetTime;
        hzPreLight = length(locsPreLight)/lightDur;
        hzDuringLight = length(locsDuringLight)/lightDur;
        hzPostLight = length(locsPostLight)/lightDur;
        isiBaseline = diff(locsBaseline);
        isiCv = std(isiBaseline)/mean(isiBaseline);
        isiMean = mean(isiBaseline);
        isiStd = std(isiBaseline);       
    end
    
    % storing sweep by sweep data in arrays for easy mean & std calculations later
    hzBaselineBySweep = [hzBaselineBySweep; hzBaseline];
    hzPreLightBySweep = [hzPreLightBySweep; hzPreLight];
    hzDuringLightBySweep = [hzDuringLightBySweep; hzDuringLight];
    hzPostLightBySweep = [hzPostLightBySweep; hzPostLight];
    allBaselineIsi = [allBaselineIsi; isiBaseline];
    isiCvBySweep = [isiCvBySweep, isiCv];
    isiMeanBySweep = [isiMeanBySweep, isiMean];
    isiStdBySweep = [isiStdBySweep, isiStd];       
    
    % create list of sweepNumber with the same size as list of timestamps 
    % to organize raster plot 
    sweepNumberArray = sweepNumber.* ones(length(locs),1);
    
    % storing sweep by sweep data in a cell array
    % to export later
    % fyi access sweep 1 data using cell2mat(tsBySweep(1))
    tsBySweep = [tsBySweep, locs];
    isiBySweep = [isiBySweep, isiBaseline];
    sweepNumberArrayBySweep = [sweepNumberArrayBySweep, sweepNumberArray];
    %----------------------------------------------------------------
    
    % checking if cell is irregular (ISI CV > 0.2)
    % ASSUMPTION ALERT, MIGHT NEED UPDATING
    if isiCv > 0.2
        isIrregular = 1;        
    elseif isiCv <= 0.2
        isIrregular = 0;        
    else
        isIrregular = NaN;        
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
    % Storing a subset of sweep by sweep data that will be exported
    data = [data; ...
        mouseNumber, ...
        experimentDate, ...
        sweepNumber, ...
        isDiscarded, ...
        gridColumns, ...
        gridRows, ...
        sweepsPerSquare, ...
        discardedSweepsFromEnd, ...
        peaksOrValleysAsNum, ...
        highpassThreshold, ...
        lowpassThreshold, ...
        minPeakHeight, ...        
        minPeakDistance, ...    
        lightExtensionFactor, ...
        lightChannel, ...
        singleLightPulse, ...
        stimDur, ... 
        stimFreq, ...
        lightDur, ...
        hzBaseline, ...
        isiMean, ...
        isiStd, ... 
        isiCv, ...
        isIrregular, ...
        length(locsPreLight), ...
        length(locsDuringLight), ...
        length(locsPostLight), ...
        hzPreLight, ...
        hzDuringLight, ...
        hzPostLight, ...
        ];   
    %----------------------------------------------------------------
    
end


%% CELL ANALYSIS - firing (all sweeps, irrespective of square) ===============================================

% Mean and Std for pre-light baseline firing rate
hzPreLightMean = mean(hzPreLightBySweep, 'omitnan');
hzPreLightStd = std(hzPreLightBySweep, 'omitnan');

% counting APs accross ALL SWEEPS
edges = [0:30];
[N, ~] = histcounts(allTimeStamps,edges);
firingHz = N/length(allSweeps);

% is firing modulated by light?
% if cell is inhibited, lightEffect = -1
% if cell is excited, lightEffect = 1
% if cell is indifferent, lightEffect = 0
lightEffectBySweep = [];
sdFromPreLightHzBySweep = [];
for sweepNumber = [1:length(allSweeps)]
    
    % calculate change in firing as standard deviations from pre-light
    % baseline mean
    sdFromPreLightHz = (hzDuringLightBySweep(sweepNumber) - hzPreLightMean) / hzPreLightStd;
    
    % store sweep-by-sweep change in firing
    sdFromPreLightHzBySweep = [sdFromPreLightHzBySweep; sdFromPreLightHz];
    
    % check if this is a discarded sweep. If so, assign NaN to stored
    % variables
    if isDiscardedBySweep(sweepNumber) == 1
        lightEffectBySweep = [lightEffectBySweep; NaN];
       
    % if this is not a discarded sweep, assign a "boolean" value to the
    % variable lightEffect. -1 means that the cell was supressed. +1 means
    % that the cell was excited. 0 means that the cell did not change its
    % firing rate
    elseif sdFromPreLightHz < -2
        lightEffectBySweep = [lightEffectBySweep; -1];
        
    elseif sdFromPreLightHz > 2
        lightEffectBySweep = [lightEffectBySweep; +1];
        
    else 
        lightEffectBySweep = [lightEffectBySweep; 0];
        
    end
end

% add lightEffect and sdFromPreLightHz as the last columns of the sweep by sweep data
data = [data, lightEffectBySweep, sdFromPreLightHzBySweep];

% Check if cell is irregular
isIrregularCell = median(isIrregularBySweep, 'omitnan');


%% CELL ANALYSIS - AP shape (all sweeps, irrespective of square)==============================================

% calculate average AP shape/trace
avgAP = mean(ySubsetAll);
%--------------------------------------------------------------------------

% create x axis for plotting AP shape
xSubset = 1000*linspace(0,(preAPinSeconds + postAPinSeconds), length(ySubsetForAPshape));
%--------------------------------------------------------------------------

% find AP peak and valley
avgAPpeakInDataPoints = find(avgAP==max(avgAP));
avgAPvalleyInDataPoints = find(avgAP==min(avgAP));
avgAPpeakInMilliSeconds = xSubset(avgAPpeakInDataPoints);
avgAPvalleyInMilliSeconds = xSubset(avgAPvalleyInDataPoints);
%--------------------------------------------------------------------------

% calculating derivatives and creating xAxis to plot derivatives. 
% for each derivative, the xAxis length decreases by 1 point.
xForDy = xSubset;
xForDy(end) = [];   % remove last point from xAxis
dy = diff(avgAP)./diff(xSubset);
ddy = diff(dy)./diff(xForDy);
xForDdy = xForDy;
xForDdy(end) = [];  % remove last point from xAxis
%--------------------------------------------------------------------------

% Find AP ONset based on 2nd derivative - ddy == min of first valley
[pks1,locs1,w,p] = findpeaks(-ddy,xForDdy,'MinPeakHeight',ddyValleyThreshold);
ddyBasedOnsetInMilliSeconds = locs1(1);       
ddyBasedOnsetInDataPoints = round(ddyBasedOnsetInMilliSeconds*(samplingFrequency/1000));
%--------------------------------------------------------------------------

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
%--------------------------------------------------------------------------

% Calculate AP width and duration based on multiple criteria
halfWidth = avgAPpeakInMilliSeconds - avgAPvalleyInMilliSeconds;
biphasicDuration = avgAPpeakInMilliSeconds - ddyBasedOnsetInMilliSeconds;
% I am commenting this duration calculation out because the
% ddyBasedOffsetInMilliSeconds is often incorrect and it takes one extra
% uneccessary user input.
% totalDurationDdyBased = ddyBasedOffsetInMilliSeconds - ddyBasedOnsetInMilliSeconds;
totalDurationAvgBased = avgAPoffsetInMilliSeconds - ddyBasedOnsetInMilliSeconds;
%--------------------------------------------------------------------------

% Check if AP duration is consistent with DA cell
% DA cells have total duration > 2 ms
% To be conservative, I'm taking the average between my two duration
% metrics: if (totalDurationDdyBased + totalDurationAvgBased)/2 > 2
% UPDATE Aug 30 2022: don't take the average. Just ignore totalDurationDdyBased.
if totalDurationAvgBased > 2
    isDA = 1;
else
    isDA = 0;
end
%--------------------------------------------------------------------------

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
    nAPtotal, ...
    ddyBasedOnsetInMilliSeconds, ...
    avgAPvalleyInMilliSeconds, ...
    avgAPpeakInMilliSeconds, ...
    avgAPoffsetInMilliSeconds, ...
    halfWidth, ...
    biphasicDuration, ...
    totalDurationAvgBased, ...
    isDA];


%% CELL ANALYSIS - POLYGON SPECIFIC CODE - data for each square ==============================================
% ALERT: I will have to adjust this code to account for removed sweeps!

% Store cell-specific data
dataCell = [mouseNumber, ...
    experimentDate, ...
    firstSweepNumber, ...
    lastSweepNumber, ...
    length(allSweeps) - length(discardedSweeps), ...
    gridColumns, ...
    gridRows, ...
    peaksOrValleysAsNum, ...
    highpassThreshold, ...
    lowpassThreshold, ...
    minPeakHeight, ...        
    minPeakDistance, ...    
    lightExtensionFactor, ...
    lightChannel, ...
    singleLightPulse, ...
    heatmapMin, ...
    heatmapMax, ...
    stimDur, ... 
    stimFreq, ...
    lightDur, ...    
    halfWidth, ...
    totalDurationAvgBased, ...    
    mean(hzBaselineBySweep, 'omitnan'), ...
    std(hzBaselineBySweep, 'omitnan'), ...
    std(hzBaselineBySweep, 'omitnan')/mean(hzBaselineBySweep, 'omitnan'), ...
    mean(allBaselineIsi, 'omitnan'), ...
    std(allBaselineIsi, 'omitnan'), ...
    std(allBaselineIsi, 'omitnan')/mean(allBaselineIsi, 'omitnan'), ...
    isDA, ...
    isIrregularCell];

% Assigning sweeps to each o-stim area (square)
% Each row in sweepIDperSquare correspond to a square
% Each column in sweepIDperSquare corresponds to a repeat sweep in that
% square.
totalSquares = gridColumns * gridRows;
relativeSweepNumberPerSquare = []; % from 1 to totalSquares
absoluteSweepNumberPerSquare = []; % from 1st sweep to last sweep number
dataSquare = [];
dataCellSingleRow = dataCell;

% to properly plot histograms later, I need to keep track of sweeps per
% square taking into consideration any discarded sweeps. I will call this
% variable sweepsPerSquarePerSquare
sweepsPerSquarePerSquareAll = [];

% each row is a square
for row=[1:totalSquares]
    
    % create/erase arrays that will be filled for each square
    dataSubsetForMeanAndSD = [];
    dataSubsetForMode = [];
    dataSubsetForMedian = [];
    
    % each column is a sweep
    for column=[1:sweepsPerSquare]    
               
        % store relative sweepID (1 to total # of sweeps)
        relativeSweepNumber = row + totalSquares * (column - 1);
        relativeSweepNumberPerSquare(row,column) = relativeSweepNumber;
        
        % store absolute sweepID (firstSweepNumber to lastSweepNumber)
        absoluteSweepNumber = relativeSweepNumber + firstSweepNumber - 1;                
        absoluteSweepNumberPerSquare(row,column) = absoluteSweepNumber;
        
        % ignore data from discarded sweeps
        if ismember(absoluteSweepNumber, discardedSweeps)
            disp('heyo I blocked trash sweep from polygon data');
        
        % save data from other sweeps    
        else
            dataSubsetForMeanAndSD = [dataSubsetForMeanAndSD; ...
                                        hzPreLightBySweep(relativeSweepNumber), ...
                                        hzDuringLightBySweep(relativeSweepNumber), ...
                                        hzPostLightBySweep(relativeSweepNumber)];
                                    
            dataSubsetForMode = [dataSubsetForMode; ...
                                        lightEffectBySweep(relativeSweepNumber)];
            
            dataSubsetForMedian = [dataSubsetForMedian; ...
                                        sdFromPreLightHzBySweep(relativeSweepNumber)];
        
        end 
        
    end
    
    
    % store number of sweeps analyzed for each square
    sweepsPerSquarePerSquare = size(dataSubsetForMeanAndSD,1);
    sweepsPerSquarePerSquareAll = [sweepsPerSquarePerSquareAll; sweepsPerSquarePerSquare];
    
    % store data for each square
    dataSquare = [dataSquare; ...
                    sweepsPerSquarePerSquare, ...
                    mean(dataSubsetForMeanAndSD), ...
                    std(dataSubsetForMeanAndSD), ...
                    mode(dataSubsetForMode), ...
                    median(dataSubsetForMedian)];
    
    % store square-by-square data together with cell data
    dataCellSingleRow = [dataCellSingleRow, ...
                    sweepsPerSquarePerSquare, ...
                    mean(dataSubsetForMeanAndSD), ...
                    std(dataSubsetForMeanAndSD), ...
                    mode(dataSubsetForMode), ...
                    median(dataSubsetForMedian)];   
end

% concatenate dataSquare with absoluteSweepNumberPerSquare and a column
% vector with square ID numbers
dataSquareWithSweeps = [[1:totalSquares]', absoluteSweepNumberPerSquare, dataSquare];

% concatenate dataSquareWithSweeps with dataCell
% use repmat to copy dataCell into multiple rows - as many rows as there
% are squares
dataCellMultipleRows = [repmat(dataCell,totalSquares,1), dataSquareWithSweeps];



%% PLOT 1 - cropped cell image ===============================================================================
% make sure you're getting the image taken with zoom = 1
% concatenate strings from user input to get full path to figure file
cellImageFileDir = strcat(cellImageDir,'\',cellImageFileName);

% import image
cellImage = imread(cellImageFileDir);

% crop image according to the polygon's mirror calibration
% I verified this in PPT lol
% the original image is 1376 x 1024 pixels
% ASSUMPTION ALERT: the calibration of the polygon will not change over time
% I need to crop the top 100 pixels and the bottom 51 pixels
% imcrop determines the rectangle to keep in the following form: [XMIN YMIN WIDTH HEIGHT]
% note that y increases from top to bottom, so ymin should match my
% required "top crop".
% I do not need to crop along the x axis, so xmin = 1 and width = 1376
% the height is 1024 - 100 - 51 = 873
croppedImage = imcrop(cellImage, [1,100,1376,873]);
figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - cell image'));
imshow(croppedImage, 'Border', 'tight');

% get inner figure size and store half of those values
pos = get(gcf, 'InnerPosition'); %// gives x left, y bottom, width, height
innerWidth = pos(3)/2;
innerHeight = pos(4)/2;

% get outer figure size and store half of those values
pos = get(gcf, 'OuterPosition'); %// gives x left, y bottom, width, height
outerWidth = pos(3)/2;
outerHeight = pos(4)/2;

% set figure size to the stored values
set(gcf,'InnerPosition',[0 maxHeight-innerHeight innerWidth innerHeight]);



%% PLOT 2 - Polygon Heatmap 1 ================================================================================

% re-organize data as a grid for heatmap
% the " .' " at the end makes sure that the sweeps are placed in the
% correct place in the grid, given an ordered polygon design - aka in each
% sweep the light moves one square to the right, or to the first square in
% the next row when it reached the end of the current row.
dataForHeatmap = reshape(dataSquare(:,end),gridColumns,[]).';

% plot heatmap as an actual heatmap
figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - heatmap 1'));

% set the ColorLimits for consistent color-coding accross cells/mice
% set the ColorMap to a flipped map with "flipud", so that more negative
% values will be "hotter" and values closer to zero will be "colder".
% remove labels from each square with CellLabelColor = 'none'.
h = heatmap(dataForHeatmap, 'ColorLimits', [heatmapMin,heatmapMax], 'ColorMap', flipud(parula), 'CellLabelColor', 'none');

% removing labels from x axis
cdl = h.XDisplayLabels;                                    % Current Display Labels
h.XDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels

% removing labels from y axis
cdl = h.YDisplayLabels;                                    % Current Display Labels
h.YDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels

% re-size 
set(gcf,'InnerPosition',[2*innerWidth maxHeight-innerHeight innerWidth innerHeight]);



%% PLOT 3 - Polygon Heatmap 2 ================================================================================
% Made to be overlayed on top of cell image from rig

% resize the heatmap matrix to match the size of the cropped image from the rig
% the original heatmap matrix will only be as big as the number of squares
% used in the polygon design. But the image from the rig will be a lot
% bigger (the size of the matrix is the size of image in pixels)
resizedHeatmap = imresize(dataForHeatmap, [size(croppedImage,1) size(croppedImage,2)], 'nearest');

% make heatmap without the heatmap function
figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - heatmap 2'));
imshow(resizedHeatmap,'Colormap',flipud(parula),'DisplayRange', [heatmapMin,heatmapMax], 'Border', 'tight');
set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
c = colorbar;
c.Label.String = 'SDs from baseline firing rate';



%% PLOT 4 - Make a blended image with transparency ===========================================================
% White where there was no inhibition of firing
% Transparent where there was inhibition of firing beyond significance
% criteria (firing decreased by -2 SDs from mean baseline firing rate)
% Note that actual criteria is determined by user input "heatmapMin"
% If heatmapMin = -2, then the criteria is "firing decreased by -2 SDs from
% mean baseline firing rate"

% pilot code that is now obsolete:
% figure('Name', '1')
% background = im2double(croppedImage);
% blendedImage = 0.5*resizedHeatmap + background;
% blendedImage = im2uint8(blendedImage);
% imshow(blendedImage)

% creating a transparency mask
% trasform all negative values into positive ones: -resizedHeatmap
transparencyMask = -resizedHeatmap;
% default values are heatmapMax=0 and heatmapMin=-2
% so default results are transparencyMaskMin=0 and transparencyMaskMax=2
transparencyMaskMin = -heatmapMax;
transparencyMaskMax = -heatmapMin;
% anything between 0 and 2 --> normalized between 1 and 0 --> transparency gradient
transparencyMask(transparencyMask>transparencyMaskMin & transparencyMask<transparencyMaskMax) = 1 - transparencyMask(transparencyMask>transparencyMaskMin & transparencyMask<transparencyMaskMax)/transparencyMaskMax;
% anything below 0 --> 1 --> min transparency --> white
transparencyMask(transparencyMask<=transparencyMaskMin) = 1;
% anything above 2 --> 0 --> max transparency --> transparent
transparencyMask(transparencyMask>=transparencyMaskMax) = 0;
% sanity check that transparencyMask is between 0 and 1
max(max(transparencyMask));
min(min(transparencyMask));

figure('Name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - transparent heatmap'))
background = imshow(croppedImage, 'Border', 'tight');
% make a white image the same size as the croppedImage
whiteImage = cat(3, ones(size(croppedImage)), ones(size(croppedImage)), ones(size(croppedImage)));
% fyi the placement of hold on and hold off in this code matter a lot for
% whatever reason. Any other placement leads to image overwritting instead
% of overlaying.
hold on;
% overlay white mask
whiteMask = imshow(whiteImage);
hold off;
% apply transparency mask to white mask
set(whiteMask, 'AlphaData', transparencyMask);

% set figure size to the same as the cropped cell image
% FYI this only adjusts the outer figure size, not the inner figure size...
set(gcf,'InnerPosition',[0 maxHeight-2.5*innerHeight innerWidth innerHeight]);
% set(gcf,'OuterPosition',[1 1 outerWidth outerHeight]);


%% PLOT 5 - ISI ==============================================================================================

% ISI CV pre-light across all sweeps
figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - baseline ISI counts'));

% edges: from 0-1s in 1ms steps
edges = [0:0.001:1]; 

histogram(allBaselineIsi, edges);
title([strcat(fileName, ' baseline ISI counts')],'Interpreter','none');
xlabel('ISI (s)');
ylabel('Counts per bin');
axis([0 1 0 ymaxIsiCV])
xticks([0 1]);
set(gcf,'Position',[maxWidth-400 maxHeight-400 400 400]);
yticks([0 ymaxIsiCV]);



%% PLOT 6 - ISI normalized ===================================================================================

% ISI CV pre-light across all sweeps
figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - baseline ISI prob'));

% edges: from 0-1s in 1ms steps
edges = [0:0.001:1]; 

histogram(allBaselineIsi, edges, 'Normalization', 'probability', 'FaceColor', [0 0 0],'FaceAlpha', 1, 'EdgeColor', 'none');  
title([strcat(fileName, ' baseline ISI prob')],'Interpreter','none');
xlabel('ISI (s)');
ylabel('Counts ber Bin / Total Counts');
axis([0 0.5 0 0.1])
xticks([0 0.5]);
yticks([0 0.1]);
set(gcf,'Position',[maxWidth-400 maxHeight-400-500 400 400]);



%% PLOT 7 - Raster plot (tiled) ==============================================================================
% Zoomed in: mean +- 2SD Hz is from short pre-light baseline

% create figure & name it
figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - tiled raster'));
t = tiledlayout(gridRows, gridColumns);

% plotting raster plot 
for square = [1:totalSquares]
    nexttile
    hold on;
    
    for sweep = [relativeSweepNumberPerSquare(square,:)]
        plot(cell2mat(tsBySweep(sweep)), ones(size(cell2mat(sweepNumberArrayBySweep(sweep))))+sweep, '|', 'Color', 'k')        
    end
    
    % zooming in and beautifying raster plot
    axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur relativeSweepNumberPerSquare(square,1)-20 relativeSweepNumberPerSquare(square,end)+20])
%     ylabel(strcat('Sweeps (', num2str(sweepsPerSquare), ')'));
%     xlabel('Time (s)');
    yticks([]);
    xticks(0:1:10);

    % adding light stim
    rectangle('Position', [lightOnsetTime relativeSweepNumberPerSquare(square,1)-20 lightDur relativeSweepNumberPerSquare(square,end)+100], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');

    % stop plotting things on this subplot
    hold off;
    
    % remove x labels from all plots except the last
    if square ~= gridRows * gridColumns
        xticklabels([]);
    end

    % flip the y-axis so that the first sweep is at the top and the last
    % sweep is at the bottom
    set(gca, 'YDir','reverse');
%     set(gcf, 'InnerPosition', [1 1 innerWidth innerHeight]);
end

t.TileSpacing = 'compact';
t.Padding = 'compact';
xlabel(t,'Time (s)')
ylabel(t,strcat('Sweeps (', num2str(sweepsPerSquare), ')'));

% set figure size to the same as the cropped cell image
% FYI this only adjusts the outer figure size, not the inner figure size...
% set(gcf,'OuterPosition',[1 1 outerWidth outerHeight]);
set(gcf,'InnerPosition',[innerWidth maxHeight-2.5*innerHeight innerWidth innerHeight]);



%% PLOT 8 - Firing histogram (tiled) =========================================================================
% mean +- 2SD Hz is from short pre-light baseline for all sweeps

% note that these show the AVERAGE firing rate accross sweeps
% my criteria for labeling a cell as "suppressed" does NOT use AVERAGE. It
% used MEDIAN. When I calculate "sdFromPreLightHz", I take the MEDIAN
% accross sweeps, not the average. HENCE the plots don't always agree with
% the variable sdFromPreLightHz. You can get a histogram that looks like a
% suppressed cell (the average crosses the mean+-2SD, but the
% sdFromPreLightHz is not -2SD.

% figuring out how many seconds in the sweep to create histogram edges
sweepSizeDataPts = size(xch2,1);
sweepSizeSec = sweepSizeDataPts/samplingFrequency;
edges = [0:sweepSizeSec];

% create figure & name it
figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - tiled hist 1'));
t = tiledlayout(gridRows, gridColumns);

% plotting histogram plot 
for square = [1:totalSquares]
    tsBySquare = [];
    nexttile
    hold on;
    
    for sweep = [relativeSweepNumberPerSquare(square,:)]
        tsBySquare = [tsBySquare; cell2mat(tsBySweep(sweep))]; 
    end
    
    [Nsquare, edges] = histcounts(tsBySquare,edges);
    firingHzSquare = Nsquare/sweepsPerSquarePerSquareAll(square);
    histogram('BinEdges', 0:sweepSizeSec, 'BinCounts', firingHzSquare, 'DisplayStyle', 'stairs', 'EdgeColor', 'k'); 

    % plot light stim as rectangle
    rectangle('Position', [lightOnsetTime 0 lightDur ymaxhist], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');

    % plot Hz mean as horizontal line
    yline(hzPreLightMean, '--');

    % plot +- 2 SD as rectangle around mean
    % [x y width height]
    rectangle('Position', [0 hzPreLightMean-(2*hzPreLightStd) 30 4*hzPreLightStd], 'FaceColor', [0 0 0 0.1], 'EdgeColor', 'none');

    axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur 0 ymaxhist])
    yticks([0 ymaxhist]);
    xticks(0:1:10);
    
    % remove y labels from all plots except the first
    if square ~= 1 
        yticklabels([]);
    end
    
    % remove x labels from all plots except the last
    if square ~= gridRows * gridColumns
        xticklabels([]);
    end
    
    hold off;    
end

t.TileSpacing = 'compact';
t.Padding = 'compact';
xlabel(t,'Time (s)')
ylabel(t,'Firing rate (Hz)');

% set figure size to the same as the cropped cell image
% FYI this only adjusts the outer figure size, not the inner figure size...
% set(gcf,'OuterPosition',[1 1 outerWidth outerHeight]);
set(gcf,'InnerPosition',[2*innerWidth maxHeight-2.5*innerHeight innerWidth innerHeight]);



%% PLOT 9 - AP width =========================================================================================

% Plot all APs and avg AP (not filtered, baseline subtracted)
% Ddy based ONset is marked with a blue arrow >
% Avg based PEAK is marked with a red arrow ^
% Avg based VALLEY is marked with a red arrow v
% Ddy based OFFset is marked with a blue arrow <
% Avg based OFFset is marked with a red arrow <
figure('name', strcat(fileName, " ", analysisDate, ' - firing_vs_light_polygon - AP width'));     
hold on;
    plot(xSubset, ySubsetAll,'Color', [0.75, 0.75, 0.75, 0.5], 'LineWidth', 0.2);
    plot(xSubset, avgAP,'Color','black','LineWidth',1.5); 
    plot(xSubset(avgAPpeakInDataPoints), max(avgAP),'^','color', 'r');
    plot(xSubset(avgAPvalleyInDataPoints), min(avgAP),'v','color', 'r');
%     plot(ddyBasedOffsetInMilliSeconds, avgAP(ddyBasedOffsetInDataPoints), '<', 'color', 'b');
    plot(ddyBasedOnsetInMilliSeconds, avgAP(ddyBasedOnsetInDataPoints), '>', 'color', 'b');
    plot(avgAPoffsetInMilliSeconds, avgAP(avgAPoffsetInDataPoints), '<', 'color', 'r');
    xlabel('Time (ms)');
    ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
    title([fileName ' AP width'],'Interpreter','none');
hold off;
set(gcf,'Position',[maxWidth-300 maxHeight-400 400 400]);

% Display y axis inverted to match how extracelullar spikes are most often displayed in the literature
set(gca, 'YDir','reverse');

% adding scale bar
xmaxHere = 1000*(preAPinSeconds+postAPinSeconds);
line([xmaxHere-2 xmaxHere],[min(avgAP) min(avgAP)],'Color','k')
line([xmaxHere xmaxHere],[min(avgAP) min(avgAP)+10],'Color','k')
text(xmaxHere-2, min(avgAP)+5, "2 ms")
text(xmaxHere-2, min(avgAP)+10, strcat("10 ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))



%% PLOT 10 - AP width with derivatives =======================================================================

% Plot the first (blue) and second (red) derivative of the avg
% Use this dor troubleshooting and adjusting ddyPeakThreshold and ddyValleyThreshold
figure('name', strcat(fileName, " ", analysisDate, ' - firing_vs_light_polygon - AP width ddy'));
hold on;
    plot(xSubset, avgAP,'Color','black','LineWidth',1);
    plot(xForDy, dy,'Color', 'b', 'LineWidth', 1);
    plot(xForDdy, ddy,'Color', 'r', 'LineWidth', 1);
    plot(ddyBasedOnsetInMilliSeconds, ddy(ddyBasedOnsetInDataPoints+1), '>', 'color', 'b');
%     plot(ddyBasedOffsetInMilliSeconds, ddy(ddyBasedOffsetInDataPoints), '<', 'color', 'b');
%     line([0 xSubset(end)],[ddyPeakThreshold ddyPeakThreshold], 'LineStyle', '-.');
    line([0 xSubset(end)],[-ddyValleyThreshold -ddyValleyThreshold], 'LineStyle', '--');
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title([fileName ' AP width ddy'],'Interpreter','none');
%     legend('avg AP', 'avg AP dy', 'avg AP ddy', 'ddy based ONset', 'ddy based OFFset', 'ddyPeakThreshold', 'ddyValleyThreshold', 'Location', 'northeast');
    legend('avg AP', 'avg AP dy', 'avg AP ddy', 'ddy based ONset', 'ddyValleyThreshold', 'Location', 'northeast');
hold off;
set(gcf,'Position',[maxWidth-300 maxHeight-400-500 400 400]);



%% EXPORTING XLS files ==========================================

% store key sweep-by-sweep data
filename = strcat(fileName, '_', analysisDate, " - firing_vs_light_polygon - sweep_by_sweep");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(data);
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
    {'mouse', ...
    'date', ...
    'sweep', ...
    'isDiscarded(1)orNot(0)', ...
    'gridColumns', ...
    'gridRows', ...
    'sweepsPerSquare', ...    
    'discardedSweepsFromEnd', ...
    'peaks(1)OrValleys(-1)', ...
    'highpassThreshold', ...
    'lowpassThreshold', ...
    'minPeakHeight', ...        
    'minPeakDistance', ...    
    'lightExtensionFactor', ...
    'lightChannel', ...
    'singleLightPulse(1)orTrain(0)', ...
    'lightPulseDur(s)', ... 
    'lightStimFreq(Hz)', ...
    'lightDur(s)', ...
    'baselineHz', ...
    'baselineIsiMean(s)', ...
    'baselineIsiSD(s)', ...
    'baselineIsiCV', ...
    'irregular', ...
    'preLightAPs', ...
    'duringLightAPs', ...
    'postLightAPs', ...
    'preLightHz', ...
    'duringLightHz', ...
    'postLightHz', ...
    'lightEffect', ...
    'sdFromPreLightHz'});
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the sweep_by_sweep xls file')
%--------------------------------------------------------------------------


% store key square-by-square data - each row is a square
filename = strcat(fileName, '_', analysisDate, " - firing_vs_light_polygon - square_by_square");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(dataCellMultipleRows);
% create cell array for variable naming
headings = cell(1, sweepsPerSquare);
for sweep = [1:sweepsPerSquare]
    headings(sweep) = {strcat('sweep', num2str(sweep))};
end
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
    {'mouse', ...
    'date', ...
    'firstSweep', ...
    'lastSweep', ...
    'nSweeps', ...
    'gridColumns', ...
    'gridRows', ...
    'peaks(1)OrValleys(-1)', ...
    'highpassThreshold', ...
    'lowpassThreshold', ...
    'minPeakHeight', ...        
    'minPeakDistance', ...
    'lightExtensionFactor', ...
    'lightChannel', ...
    'singleLightPulse(1)orTrain(0)', ...
    'heatmapMin', ...
    'heatmapMax', ...
    'lightPulseDur(s)', ... 
    'lightStimFreq(Hz)', ...
    'lightDur(s)', ...   
    'APhalfWidth(ms)', ...
    'APduration(ms)', ...
    'baselineHzMean', ...
    'baselineHzSD', ...
    'baselineHzCV', ...
    'baselineISImean', ...
    'baselineISIsd', ...
    'baselineISIcv', ...
    'isDA(0,1)', ...
    'isIrregular(0,1)', ...
    'square', ...
    headings{:}, ...
    'sweepsPerSquare', ...
    'preLight(Hz)AVG', ...
    'duringLight(Hz)AVG', ...
    'postLight(Hz)AVG', ...
    'preLight(Hz)STD', ...
    'duringLight(Hz)STD', ...
    'postLight(Hz)STD', ...
    'lightEffect(-1,0,1)', ...
    'medianSdFromPreLightHz'});
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the square_by_square xls file')
%--------------------------------------------------------------------------


% store all timestamps from all sweeps in a single column
filename = strcat(fileName, '_', analysisDate, " - firing_vs_light_polygon - all_AP_timestamps");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
writematrix(allTimeStamps, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the all_AP_timestamps xls file')
%--------------------------------------------------------------------------


% store AP shape data in a single row
filename = strcat(fileName, '_', analysisDate, " - firing_vs_light_polygon - AP_shape");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(dataAPshape);
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
    {'mouse', ...
    'date', ...
    'firstSweep', ...
    'lastSweep', ...
    'discardedSweepsFromEnd', ...
    'peaks(1)OrValleys(-1)', ...
    'highpassThreshold', ...
    'lowpassThreshold', ...
    'minPeakHeight', ...        
    'minPeakDistance', ...
    'preAP(s)', ...
    'postAP(s)', ...
    'baselinePreAPdur(s)', ...
    'ddyValleyThreshold', ...
    'nAPtotal', ...
    'ddyBasedOnset(ms)', ...
    'APvalley(ms)', ...
    'APpeak(ms)', ...
    'avgBasedOffset(ms)', ...
    'halfWidth(ms)', ...
    'biphasicDuration(ms)', ...
    'totalDurationAvgBased(ms)', ...
    'isDA'});
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the AP_shape xls file')
%--------------------------------------------------------------------------


% store cell data with square-by-square info in a single row
filename = strcat(fileName, '_', analysisDate, " - firing_vs_light_polygon - square_by_square cell");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(dataCellSingleRow);
dataInCellFormat = [dataInCellFormat, cellImageFileName];

% create cell array with strings for naming square-by-square data
squareVariables = cell(1,size(dataSquare,1)*size(dataSquare,2));
for square = [1:totalSquares]
    startingIndex = square*9 - 8;
    squareVariables(startingIndex) = {strcat('sq', num2str(square), '-sweepsPerSquare')};
    squareVariables(startingIndex+1) = {strcat('sq', num2str(square), '-preLight(Hz)AVG')};
    squareVariables(startingIndex+2) = {strcat('sq', num2str(square), '-duringLight(Hz)AVG')};
    squareVariables(startingIndex+3) = {strcat('sq', num2str(square), '-postLight(Hz)AVG')};
    squareVariables(startingIndex+4) = {strcat('sq', num2str(square), '-preLight(Hz)STD')};
    squareVariables(startingIndex+5) = {strcat('sq', num2str(square), '-duringLight(Hz)STD')};
    squareVariables(startingIndex+6) = {strcat('sq', num2str(square), '-postLight(Hz)STD')};
    squareVariables(startingIndex+7) = {strcat('sq', num2str(square), '-lightEffect(-1,0,1)')};
    squareVariables(startingIndex+8) = {strcat('sq', num2str(square), '-medianSdFromPreLightHz')};
end

% label stored data
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
    {'mouse', ...
    'date', ...
    'firstSweep', ...
    'lastSweep', ...
    'nSweeps', ...
    'gridColumns', ...
    'gridRows', ...
    'peaks(1)OrValleys(-1)', ...
    'highpassThreshold', ...
    'lowpassThreshold', ...
    'minPeakHeight', ...        
    'minPeakDistance', ...
    'lightExtensionFactor', ...
    'lightChannel', ...
    'singleLightPulse(1)orTrain(0)', ...
    'heatmapMin', ...
    'heatmapMax', ...
    'lightPulseDur(s)', ... 
    'lightStimFreq(Hz)', ...
    'lightDur(s)', ...   
    'APhalfWidth(ms)', ...
    'APduration(ms)', ...
    'baselineHzMean', ...
    'baselineHzSD', ...
    'baselineHzCV', ...
    'baselineISImean', ...
    'baselineISIsd', ...
    'baselineISIcv', ...
    'isDA(0,1)', ...
    'isIrregular(0,1)', ...
    squareVariables{:}, ...
    'cellImageFileName'});
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the square_by_square cell xls file')



%% Obsolete code

%% PLOT 7.1 - Raster plot (subplots) ================================================================================
% Obsolete code - I replaced subplots by tiled plots
% Zoomed in: mean +- 2SD Hz is from short pre-light baseline

% % create figure & name it
% figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - raster'));
% 
% % plotting raster plot 
% for square = [1:totalSquares]
%     subplot(gridRows,gridColumns,square)
%     hold on;
%     
%     for sweep = [sweepOrderPerSquare(square,:)]
%         plot(cell2mat(tsBySweep(sweep)), ones(size(cell2mat(sweepNumberArrayBySweep(sweep))))+sweep, '|', 'Color', 'k')        
%     end
%     
%     % zooming in and beautifying raster plot
%     axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur sweepOrderPerSquare(square,1)-20 sweepOrderPerSquare(square,end)+20])
%     ylabel(strcat('Sweeps (', num2str(sweepsPerSquare), ')'));
%     xlabel('Time (s)');
%     yticks([]);
%     xticks(0:1:10);
%     xticklabels([]);
% 
%     % adding light stim
%     rectangle('Position', [lightOnsetTime sweepOrderPerSquare(square,1)-20 lightDur sweepOrderPerSquare(square,end)+100], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');
% 
%     % stop plotting things on this subplot
%     hold off;
% 
%     % flip the y-axis so that the first sweep is at the top and the last
%     % sweep is at the bottom
%     set(gca, 'YDir','reverse');
%     set(gcf, 'InnerPosition', [0 maxHeight-innerHeight-400 innerWidth innerHeight]);
% end


%% PLOT 8.1 - Firing histogram (subplots) ======================================
% Obsolete code - I replaced subplots by tiled plots
% mean +- 2SD Hz is from short pre-light baseline for all sweeps

% % note that these show the AVERAGE firing rate accross sweeps
% % my criteria for labeling a cell as "suppressed" does NOT use AVERAGE. It
% % used MEDIAN. When I calculate "sdFromPreLightHz", I take the MEDIAN
% % accross sweeps, not the average. HENCE the plots don't always agree with
% % the variable sdFromPreLightHz. You can get a histogram that looks like a
% % suppressed cell (the average crosses the mean+-2SD, but the
% % sdFromPreLightHz is not -2SD.
% 
% % figuring out how many seconds in the sweep to create histogram edges
% sweepSizeDataPts = size(xch2,1);
% sweepSizeSec = sweepSizeDataPts/samplingFrequency;
% edges = [0:sweepSizeSec];
% 
% figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - hist 1'));
% 
% % plotting histogram plot 
% for square = [1:totalSquares]
%     tsBySquare = [];
%     subplot(gridRows,gridColumns,square)
% 
%     for sweep = [sweepOrderPerSquare(square,:)]
%         tsBySquare = [tsBySquare; cell2mat(tsBySweep(sweep))]; 
%     end
% 
%     [Nsquare, edges] = histcounts(tsBySquare,edges);
%     firingHzSquare = Nsquare/sweepsPerSquare;
%     histogram('BinEdges', 0:sweepSizeSec, 'BinCounts', firingHzSquare, 'DisplayStyle', 'stairs', 'EdgeColor', 'k'); 
% 
%     % plot light stim as rectangle
%     rectangle('Position', [lightOnsetTime 0 lightDur ymaxhist], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');
% 
%     % plot Hz mean as horizontal line
%     yline(hzPreLightMean, '--');
% 
%     % plot +- 2 SD as rectangle around mean
%     % [x y width height]
%     rectangle('Position', [0 hzPreLightMean-(2*hzPreLightStd) 30 4*hzPreLightStd], 'FaceColor', [0 0 0 0.1], 'EdgeColor', 'none');
% 
%     xlabel('Time (s)');
%     ylabel('Firing rate (Hz)');
%     axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur 0 ymaxhist])
%     set(gcf, 'InnerPosition', [innerWidth maxHeight-innerHeight innerWidth innerHeight]);
%     yticks([0 ymaxhist]);
%     hold off;     
% end


%% PLOT 9.1 - Firing histogram 2 (subplots) =============================
% Obsolete code - I replaced subplots by tiled plots
% Obsolete code - I prefer to compare the square-by-square data to the
% cell's average firing rate and not the square-specific firing rate
% mean +- 2SD Hz is from short pre-light baseline for specific squares

% % figuring out how many seconds in the sweep to create histogram edges
% sweepSizeDataPts = size(xch2,1);
% sweepSizeSec = sweepSizeDataPts/samplingFrequency;
% edges = [0:sweepSizeSec];
% 
% figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - hist 2'));
% 
% % plotting histogram plot 
% for square = [1:totalSquares]
%     tsBySquare = [];
%     subplot(gridRows,gridColumns,square)
% 
%     for sweep = [sweepOrderPerSquare(square,:)]
%         tsBySquare = [tsBySquare; cell2mat(tsBySweep(sweep))]; 
%     end
% 
%     [Nsquare, edges] = histcounts(tsBySquare,edges);
%     firingHzSquare = Nsquare/sweepsPerSquare;
%     histogram('BinEdges', 0:sweepSizeSec, 'BinCounts', firingHzSquare, 'DisplayStyle', 'stairs', 'EdgeColor', 'k'); 
% 
%     % plot light stim as rectangle
%     rectangle('Position', [lightOnsetTime 0 lightDur ymaxhist], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');
% 
%     % plot Hz mean as horizontal line
%     hzPreLightMeanSquare = dataSquare(square,1);
%     yline(hzPreLightMeanSquare, '--');
% 
%     % plot +- 2 SD as rectangle around mean
%     % [x y width height]
%     hzPreLightStdSquare = dataSquare(square,4);
%     rectangle('Position', [0 hzPreLightMeanSquare-(2*hzPreLightStdSquare) 30 4*hzPreLightStdSquare], 'FaceColor', [0 0 0 0.1], 'EdgeColor', 'none');
% 
%     xlabel('Time (s)');
%     ylabel('Firing rate (Hz)');
%     axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur 0 ymaxhist])
%     set(gcf, 'InnerPosition', [innerWidth maxHeight-innerHeight-100 innerWidth innerHeight]);
%     yticks([0 ymaxhist]);
%     hold off;     
% end


%% PLOT 9.2 - Firing histogram 2 (tiled) ================================================================================
% Obsolete code - I prefer to compare the square-by-square data to the
% cell's average firing rate and not the square-specific firing rate
% mean +- 2SD Hz is from short pre-light baseline for specific squares

% % figuring out how many seconds in the sweep to create histogram edges
% sweepSizeDataPts = size(xch2,1);
% sweepSizeSec = sweepSizeDataPts/samplingFrequency;
% edges = [0:sweepSizeSec];
% 
% % create figure & name it
% figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - tiled hist 2'));
% t = tiledlayout(gridRows, gridColumns);
% 
% % plotting histogram plot 
% for square = [1:totalSquares]
%     tsBySquare = [];
%     nexttile
%     hold on;
%     
%     for sweep = [sweepOrderPerSquare(square,:)]
%         tsBySquare = [tsBySquare; cell2mat(tsBySweep(sweep))]; 
%     end
%     
%     [Nsquare, edges] = histcounts(tsBySquare,edges);
%     firingHzSquare = Nsquare/sweepsPerSquare;
%     histogram('BinEdges', 0:sweepSizeSec, 'BinCounts', firingHzSquare, 'DisplayStyle', 'stairs', 'EdgeColor', 'k'); 
% 
%     % plot light stim as rectangle
%     rectangle('Position', [lightOnsetTime 0 lightDur ymaxhist], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');
% 
%     % plot Hz mean as horizontal line
%     hzPreLightMeanSquare = dataSquare(square,1);
%     yline(hzPreLightMeanSquare, '--');
% 
%     % plot +- 2 SD as rectangle around mean
%     % [x y width height]
%     hzPreLightStdSquare = dataSquare(square,4);
%     rectangle('Position', [0 hzPreLightMeanSquare-(2*hzPreLightStdSquare) 30 4*hzPreLightStdSquare], 'FaceColor', [0 0 0 0.1], 'EdgeColor', 'none');
% 
%     axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur 0 ymaxhist])
%     yticks([0 ymaxhist]);
%     xticks(0:1:10);
%     
%     % remove y labels from all plots except the first
%     if square ~= 1 
%         yticklabels([]);
%     end
%     
%     % remove x labels from all plots except the last
%     if square ~= gridRows * gridColumns
%         xticklabels([]);
%     end
%     
%     hold off;    
% end
% 
% t.TileSpacing = 'compact';
% t.Padding = 'compact';
% xlabel(t,'Time (s)')
% ylabel(t,'Firing rate (Hz)');
% 
% % set figure size to the same as the cropped cell image
% set(gcf,'InnerPosition',[0 maxHeight-innerHeight-500 innerWidth innerHeight]);


%% PLOT - AP width with derivatives normalized ==================

% % Plot the first (blue) and second (red) derivative of the avg - Normalized
% figure('name', strcat(fileName,' y dy ddy Normalized'));
% hold on;
%     plot(xSubset, avgAP/max(avgAP),'Color','black','LineWidth',1);
%     plot(xForDy, dy/max(dy),'Color', 'b', 'LineWidth', 1);
%     plot(xForDdy, ddy/max(ddy),'Color', 'r', 'LineWidth', 1);
%     plot(ddyAfterLastPeakOrValleyInMilliSeconds, ddy(ddyAfterLastPeakOrValleyInDataPoints)/max(ddy), 'o', 'color', 'k');
%     plot(ddyBasedOffsetInMilliSeconds, ddy(ddyBasedOffsetInDataPoints)/max(ddy), '<', 'color', 'b');
%     plot(ddyBasedOnsetInMilliSeconds, ddy(ddyBasedOnsetInDataPoints)/max(ddy), '>', 'color', 'b');
%     xlabel('Time (ms)');
%     ylabel('Normalized Amplitude (au)');
%     title([fileName ' Normalized y dy ddy'],'Interpreter','none');
%     legend('avg AP', 'avg AP dy', 'avg AP ddy', 'Location', 'northwest');
% hold off;
% movegui('northeast');


%% PLOT - niceplot of first sweep ===============================

% % plotting niceplot of first sweep    
% figure('name', strcat(fileName, " ", analysisDate, ' - firing_vs_light - niceplot 1st sweep'));
% plot(x,yFilteredAll(:,1),'k','LineWidth',1);
% axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur -ymax ymax]);
% title([fileName ' - firing_vs_light - niceplot 1st sweep'],'Interpreter','none');
% set(gca,'Visible','off');
% set(gcf,'Position',[1400 50 500 400]);
% 
% % if the light stim is a train (singleLightPulse = 0), add light
% % stim cartoon as a train.
% % adding light stim as train - faithful to light stim parameters
% % note that this code will use light stim parameters from the last sweep!
% % if light stim is not the same accross all sweeps, this will be
% % misleading!
% if singleLightPulse == 0
%     for nStim=1:length(lightPulseStart)
%         line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[ymax,ymax],'Color',[0 0.4470 0.7410],'LineWidth',10)
%     end
% end
% 
% % adding light stim as rectangle
% % note that this code will use light stim parameters from the last sweep!
% % if light stim is not the same accross all sweeps, this will be
% % misleading!
% line([lightOnsetTime,lightOnsetTime+lightDur],[ymax-10,ymax-10],'Color',[0 0.4470 0.7410],'LineWidth',10)
% 
% % adding scale bar
% xmin = lightOnsetTime-lightDur;
% xmax = lightOnsetTime+2*lightDur;
% ymin = -ymax;
% line([xmax-1 xmax],[ymax ymax],'Color','k')
% line([xmax xmax],[ymax ymax-((ymax-ymin)/10)],'Color','k')
% text(xmax-1, ymax-((ymax-ymin)/20), "1 s")
% text(xmax-1, ymax-((ymax-ymin)/10), strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
% hold off;


%% PLOT - niceplot of all sweeps ================================

% % plotting niceplot     
% figure('name', strcat(fileName, " ", analysisDate, ' - firing_vs_light - niceplot all'));
% plot(x,yFilteredAll,'Color',[0, 0, 0, 0.25]);
% axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur -ymax ymax]);
% title([fileName ' - firing_vs_light - niceplot all'],'Interpreter','none');
% set(gca,'Visible','off');
% set(gcf,'Position',[1400 550 500 400]);
% 
% % adding light stim   
% % if the light stim is a train (singleLightPulse = 0), add light
% % stim cartoon as a train.
% % note that this code will use light stim parameters from the last sweep!
% % if light stim is not the same accross all sweeps, this will be
% % misleading!
% if singleLightPulse == 0
%     for nStim=1:length(lightPulseStart)
%         line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[ymax,ymax],'Color',[0 0.4470 0.7410],'LineWidth',10)
%     end
% end
% 
% % adding light stim
% % note that this code will use light stim parameters from the last sweep!
% % if light stim is not the same accross all sweeps, this will be
% % misleading!
% line([lightOnsetTime,lightOnsetTime+lightDur],[ymax-10,ymax-10],'Color',[0 0.4470 0.7410],'LineWidth',10)
% 
% % adding scale bar
% xmin = lightOnsetTime-lightDur;
% xmax = lightOnsetTime+2*lightDur;
% ymin = -ymax;
% line([xmax-1 xmax],[ymax ymax],'Color','k')
% line([xmax xmax],[ymax ymax-((ymax-ymin)/10)],'Color','k')
% text(xmax-1, ymax-((ymax-ymin)/20), "1 s")
% text(xmax-1, ymax-((ymax-ymin)/10), strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))


%     % storing sweep by sweep data in a structure
%     % fyi access sweep 1 data using dataPerSweep.s0128ts
%     dataPerSweep.(strcat('s',num2str(sweepNumber),'ts')) = locs;
%     dataPerSweep.(strcat('s',num2str(sweepNumber),'hz')) = length(locsBaseline)/lightOnsetTime;
%     dataPerSweep.(strcat('s',num2str(sweepNumber),'isi')) = isiBaseline;


% Obsolete AP shape code 
% % Find AP OFFset based on 2nd derivative - ddy == zero after last peak
% [pks2,locs2,w,p] = findpeaks(ddy,xForDdy,'MinPeakHeight',ddyPeakThreshold);
% ddyLastPeakInMilliSeconds = locs2(end);     
% ddyLastPeakInDataPoints = round(ddyLastPeakInMilliSeconds*(samplingFrequency/1000));
% ddyLastValleyInMilliSeconds = locs1(end);     
% ddyLastValleyInDataPoints = round(ddyLastValleyInMilliSeconds*(samplingFrequency/1000));
% ddyAfterLastPeakOrValley = ddy;
% ddyBasedOffsetInDataPoints = [];
% 
% % if last peak is before last valley, look for when ddy crosses 0 after last valley
% if ddyLastPeakInMilliSeconds < ddyLastValleyInMilliSeconds
%     ddyAfterLastPeakOrValley(1:ddyLastValleyInDataPoints) = [];
%     ddyAfterLastPeakOrValleyInDataPoints = ddyLastValleyInDataPoints;
%     ddyCrossesZeroPt = find(diff(sign(ddyAfterLastPeakOrValley))>0, 1);
%     ddyBasedOffsetInDataPoints = ddyLastValleyInDataPoints + ddyCrossesZeroPt + 1;
% % if last peak is after last valley, look for when ddy crosses 0 after last peak
% else
%     ddyAfterLastPeakOrValley(1:ddyLastPeakInDataPoints) = [];
%     ddyAfterLastPeakOrValleyInDataPoints = ddyLastPeakInDataPoints;
%     ddyCrossesZeroPt = find(diff(sign(ddyAfterLastPeakOrValley))<0, 1);
%     ddyBasedOffsetInDataPoints = ddyLastPeakInDataPoints + ddyCrossesZeroPt + 1;
% end 
% 
% ddyBasedOffsetInMilliSeconds = xForDdy(ddyBasedOffsetInDataPoints);
% ddyAfterLastPeakOrValleyInMilliSeconds = xForDdy(ddyAfterLastPeakOrValleyInDataPoints);