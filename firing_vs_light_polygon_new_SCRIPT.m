%{ 
DOCUMENTATION
Created: 2022 11 04
Works? Yes
Author: PA

This function is used to analyze light-evoked changes in firing rate in
sCRACM experiments using the polygon.
This function was derived from firing_vs_light_polygon, but went through A LOT
of rewritting based on the function psc_vs_light_polygon_new --> The data
organization is very different, to better handle discarded sweeps without
losing track of which sweeps belong to which ROIs (Regions Of Interest -
subarea of the field that was illuminated). 
FYI I sometimes use the words "square" and "ROI" interchangeably.

INPUTS explained:
    - gridColumns: number of columns in the polygon ROI grid

    - gridRows: number of rows in the polygon ROI grid

    - orderedGrid: booleaen variable (0 or 1) determining whether you used
    an ordered grid or not using the polygon. In an ordered grid, ROIs are
    stimulated in order - the o-stim moves from the top left to the
    bottom right. For instance, if you have an ordered 3x3 grid, this would be the
    order of the ROIs:
     1     2     3
     4     5     6
     7     8     9

    - orderOfROIs: column representing the order of the illuminated ROIs
    when using a non-ordered grid. When piloting sCRACM, I used ordered
    grids to make on-line and off-line analysis easier. However, the gold
    standard for these experiments is to randomize the o-stim order. I
    decided to use a pseudo-random order that maximizes the distance
    between each o-stim. I called this design "spaced out". The orderOfROIs
    of the ordered 3x3 grid shown above would be:
    [1 2 3 4 5 6 7 8 9]'
    (the little apostrophe at the end transposes the row into a column)
    The orderOfROIs of the design "5x5 spaced out" is:
    [8 16 14 23 3 10 12 25 21 19 6 18 1 5 11 4 22 15 13 7 2 20 24 17 9]' 

    - discardedSweeps: specific sweeps that should NOT be analyzed due to
    artifacts. If I want to discard sweep 0024, just write 24.

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
    gridColumns = 5;
    gridRows = 5;
    orderedGrid = 0;       % 0 if NOT ordered, 1 if ordered
    orderOfROIs = [8 16 14 23 3 10 12 25 21 19 6 18 1 5 11 4 22 15 13 7 2 20 24 17 9]';     % this is the order of the design 5x5 spaced out

    discardedSweeps = [];
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
    - Width and frequency of polygon stim is the the same as LED o-stim
    (only relevant if you use the polygon channel as the light stim
    channel)

BEWARE:
    - if you add more variables into "data" for exporting, you need to
    adjust the code at multiple places. i think I finally fixed this

TO DO:
    - test with light train
%}

obj = m570.s0076;

% function firing_vs_light_polygon_new(obj)
%%  USER INPUT ==================================================

% Affects data analysis - Organizing data by o-stim grid
gridColumns = 5;
gridRows = 5;
orderedGrid = 0;       % 0 if NOT ordered, 1 if ordered
orderOfROIs = [8 16 14 23 3 10 12 25 21 19 6 18 1 5 11 4 22 15 13 7 2 20 24 17 9]';     % this is the order of the design 5x5 spaced out

% Affects data analysis - Finding APs:
discardedSweeps = [134:137];
peaksOrValleys = 'v';   
highpassThreshold = 100;
lowpassThreshold = 1500;    
minPeakHeight = 25;         
minPeakDistance = 0.05; 
lightExtensionFactor = 1;
lightChannel = 2;
singleLightPulse = 1; 

% Affects data analysis - AP shape:
preAPinSeconds = 0.005;            
postAPinSeconds = 0.01;           
preAPbaselineDurationSeconds = 0.002;
ddyValleyThreshold = 600;
  
% Affects data display: 
ymax = 100;
ymaxhist = 15;
zoomWindow = 0.25;
ymaxIsiCV = 150;
heatmapMin = -2;
heatmapMax = 0;
cellImageFileNameDIC = 's1c1_z1_dic.tif';
cellImageFileNameAlexa = 's1c1_z1_647_SUM_Stack.tif';
cellImageDir = 'E:\Priscilla - BACKUP 20200319\Ephys\2022\20221102 m570 dat nphr';

% Affects data saving:
savefileto = 'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries\2022\2022-11-03 dat nphr spaced out grid';


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

% copied from psc_vs_light_polygon_new
% ALERT: NEED TO TEST
% checking for incomplete sweeps and adding them to the list of discarded sweeps
% to avoid this error: "Index exceeds the number of array elements (0)". 
if numel(fieldnames(obj.sweeps)) < obj.header.NSweepsPerRun  
    lastCompleteSweep = firstSweepNumber + numel(fieldnames(obj.sweeps)) - 1;
    discardedSweeps = [discardedSweeps, lastCompleteSweep:lastSweepNumber];
end

% copied from psc_vs_light_polygon_new
% organizing sweeps according to total ROIs and discarding discardedSweeps
totalROIs = gridColumns * gridRows;
plannedSweepsPerROI = size(allSweeps,2)/totalROIs;
% reorganize allSweeps into an array with "totalROIs" rows and "sweepsPerROI" columns
allSweepsReorganized = reshape(allSweeps, totalROIs, plannedSweepsPerROI);
% make cell array with Sweeps/ROI info
% use a cell array so I can assign a different number of sweeps per ROI if needed
sweepsInROI = cell(totalROIs,1);
relativeSweepsInROI = cell(totalROIs,1);
sweepsPerROI = zeros(totalROIs,1);
counter = 0;
for row=1:totalROIs
    % only keep the sweeps that are not discarded
    % setdiff(a,b) keeps the sweeps that are in a, but not b
    sweeps = setdiff(allSweepsReorganized(row,:), discardedSweeps);
    sweepsInROI(row)={sweeps};

    % relativeSweepsInROI goes from 1 to max number of ANALYZED sweeps,
    % not the total number of recorded sweeps!
    % this is important for plotting later
    relativeSweeps = 1:size(sweeps,2);
    relativeSweeps = counter + relativeSweeps;
    counter = counter + size(sweeps,2);
    relativeSweepsInROI(row) = {relativeSweeps};
    
    % keep track of how many sweeps are in each ROI
    sweepsPerROI(row) = size(sweeps,2);
end

% copied from psc_vs_light_polygon_new
% re-ordering sweeps into the appropriate ROIs
% if the grid is ordered, the first sweep corresponds to the top left ROI
% and the last sweep corresponds to the bottom right ROI.
% the code in general assumes that the grid IS ordered. If it is NOT, you
% need to adjust it.
% orderedGrid = 0 if the grid is NOT ordered
if orderedGrid == 0
    sweepsInROI = sweepsInROI(orderOfROIs);   
    % if the first number in orderOfROIs is 8, matlab will move the 8th
    % item in sweepsInROI to the first row. If the 13th number of
    % orderOfROIs is 1, matlab will move the first row of
    % sweepsInROI to the 13th row.   
    % DO NOT re-order relativeSweepsInROI - I made this mistake earlier
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
dataPerROI = [];

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


%% ROI BY ROI AND SWEEP BY SWEEP ANALYSIS ===================================================================================

% get data from sweeps in file (only the subset we will analyze)
for ROI = 1:totalROIs  
    
    column = 0;
    
    % creating/clearing matrixes for each ROI
    dataSubsetForMeanAndSD = [];
    dataSubsetForMode = [];
    dataSubsetForMedian = [];
    
    for sweepNumber = cell2mat(sweepsInROI(ROI))
        
        column = column + 1;
           
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
            ROI, ...
            gridColumns, ...
            gridRows, ...
            orderedGrid, ...
            plannedSweepsPerROI, ...
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
       
        % store ROI data in a cell array
        % each ROI is stored in a row
        % each sweep within a ROI is stored in a column        
        hzPreLightInROI(ROI, column) = {hzPreLight};
        hzDuringLightInROI(ROI, column) = {hzDuringLight};
        hzPostLightInROI(ROI, column) = {hzPostLight};
                
    end    
    
    % store some stats per ROI
    meanHzPreLightPerROI(ROI) = mean(cell2mat(hzPreLightInROI(ROI,:)), 'omitnan');
    meanHzDuringLightPerROI(ROI) = mean(cell2mat(hzDuringLightInROI(ROI,:)), 'omitnan');
    meanHzPostLightPerROI(ROI) = mean(cell2mat(hzPostLightInROI(ROI,:)), 'omitnan');
    
    stdHzPreLightPerROI(ROI) = std(cell2mat(hzPreLightInROI(ROI,:)), 'omitnan');
    stdHzDuringLightPerROI(ROI) = std(cell2mat(hzDuringLightInROI(ROI,:)), 'omitnan');
    stdHzPostLightPerROI(ROI) = std(cell2mat(hzPostLightInROI(ROI,:)), 'omitnan');
    
end


%% CELL ANALYSIS - firing (all sweeps, irrespective of ROI) ===============================================

% Mean and Std for pre-light baseline firing rate
hzPreLightMean = mean(hzPreLightBySweep, 'omitnan');
hzPreLightStd = std(hzPreLightBySweep, 'omitnan');

% creating a new variable to avoid errors
totalSweepsAnalyzed = length(allSweeps) - length(discardedSweeps);

% counting APs accross ALL SWEEPS
edges = [0:30];
[N, ~] = histcounts(allTimeStamps,edges);
firingHz = N/totalSweepsAnalyzed;

% is firing modulated by light?
% if cell is inhibited, lightEffect = -1
% if cell is excited, lightEffect = 1
% if cell is indifferent, lightEffect = 0
lightEffectBySweep = [];
sdFromPreLightHzBySweep = [];
for sweepNumber = 1:totalSweepsAnalyzed
    
    % calculate change in firing as standard deviations from pre-light
    % baseline mean
    sdFromPreLightHz = (hzDuringLightBySweep(sweepNumber) - hzPreLightMean) / hzPreLightStd;
    
    % store sweep-by-sweep change in firing
    sdFromPreLightHzBySweep = [sdFromPreLightHzBySweep; sdFromPreLightHz];
         
    % assign a "boolean" value to the variable lightEffect. -1 means that
    % the cell was supressed. +1 means that the cell was excited. 0 means
    % that the cell did not change its firing rate
    if sdFromPreLightHz < -2
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
        

%% CELL ANALYSIS - AP shape (all sweeps, irrespective of ROI)==============================================

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


%% CELL ANALYSIS - POLYGON SPECIFIC CODE - data for each ROI ==============================================

% Store cell-specific data
dataCell = [mouseNumber, ...
    experimentDate, ...
    firstSweepNumber, ...
    lastSweepNumber, ...
    length(allSweeps) - length(discardedSweeps), ...
    gridColumns, ...
    gridRows, ...
    orderedGrid, ...
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

% store ROI-by-ROI info
for ROI = 1:totalROIs
    column = 0;
    
    for relativeSweepNumber = cell2mat(relativeSweepsInROI(ROI))
        column = column + 1;
        lightEffectInROI(ROI, column) = {lightEffectBySweep(relativeSweepNumber)};
        sdFromPreLightHzInROI(ROI, column) = {sdFromPreLightHzBySweep(relativeSweepNumber)};
    end
    
    lightEffectPerROI(ROI) = mode(cell2mat(lightEffectInROI(ROI, :)));
    sdFromPreLightHzPerROI(ROI) = median(cell2mat(sdFromPreLightHzInROI(ROI, :)), 'omitnan');
      
end

% store ROI-by-ROI info
dataPerROI = [meanHzPreLightPerROI', ...
    meanHzDuringLightPerROI', ...
    meanHzPostLightPerROI', ...    
    stdHzPreLightPerROI', ...
    stdHzDuringLightPerROI', ...
    stdHzPostLightPerROI', ...
    lightEffectPerROI', ...
    sdFromPreLightHzPerROI'];

% store ROI-by-ROI data in multiple rows
dataCellMultipleRows = [repmat(dataCell,totalROIs,1), ...
    [1:totalROIs]', ...
    sweepsPerROI, ...
    dataPerROI];


%% PREP for PLOTs ===================================================================================================

% make a color map from white to pink (reddish purple used in my paper)
% for coloring success rate heatmap
customColorMapPink = [linspace(1,204/255,256)' linspace(1,121/255,256)'  linspace(1,167/255,256)'];

% make color map from white to vermillion
% for coloring amplitude heatmap
customColorMapVermillion = [linspace(1,213/255,256)' linspace(1,94/255,256)'  linspace(1,0/255,256)'];

% make color map from white to blue
% for coloring onset latency heatmap
customColorMapBlue = [linspace(1,0/255,256)' linspace(1,114/255,256)'  linspace(1,178/255,256)'];

% make color map from white to sky blue
% for coloring onset latency jitter heatmap
customColorMapSky = [linspace(1,86/255,256)' linspace(1,180/255,256)'  linspace(1,233/255,256)'];


%% PLOT 1 - cropped cell image ===============================================================================
% make sure you're getting the image taken with zoom = 1
% concatenate strings from user input to get full path to figure file
cellImageFileDir = strcat(cellImageDir,'\',cellImageFileNameDIC);

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
figure('name', strcat(fileName, '_', analysisDate, ' - psc_vs_light_polygon - DIC image grid'));
hold on;
imshow(croppedImage, 'Border', 'tight');

% calculate parameters for scale bar
% ASSUMPTION ALERT: pixelsPerMicron corresponds to my usual 40x objective at 1x zoom
xmaxImg = size(croppedImage,2);    % in pixels, should be 1376
ymaxImg = size(croppedImage,1);    % in pixels, should be 874
scaleDistanceFromEdge = 50;     % in pixels
scaleBarSize = 50;              % in um
pixelsPerMicron = 873 / 222.2;  % 222.2 um in 873 pixels
scaleBarSizeInPixels = scaleBarSize * pixelsPerMicron;

% add polygon grid to figure
for row = 1:gridRows+1
    ypos = (row-1)*ymaxImg/gridRows;
    yline(ypos, 'Color', 'k', 'LineWidth', 0.5);
end
for col = 1:gridColumns+1
    xpos = (col-1)*xmaxImg/gridColumns;
    xline(xpos, 'Color', 'k', 'LineWidth', 0.5);
end

% add scale bar to figure
% line([x x], [y y])
line([scaleDistanceFromEdge scaleDistanceFromEdge+scaleBarSizeInPixels],...
    [ymaxImg-scaleDistanceFromEdge ymaxImg-scaleDistanceFromEdge],...
    'LineWidth', 2, 'Color', 'k');
% text(x, y, string)
text(scaleDistanceFromEdge,...
    ymaxImg-2*scaleDistanceFromEdge,...
    strcat(num2str(scaleBarSize), " μm"),...
    'FontSize', 10);
hold off;

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


%% PLOT 2 - cropped Alexa cell image with grid ===============================================================================
% make sure you're getting the image taken with zoom = 1
% concatenate strings from user input to get full path to figure file
cellImageFileDir = strcat(cellImageDir,'\',cellImageFileNameAlexa);

% import image
cellImage = imread(cellImageFileDir);

% normalize image to max intensity value
% this is important for dealing with summed z-stacks (instead of max
% intensity z-stacks)
if max(max(cellImage))/256 > 1
    adjustmentFactor = (max(max(cellImage))/256) * 256;
    cellImage = cellImage/adjustmentFactor;
end
    
% invert alexa image so that black = white
% I wanted to do this anyway, but MATLAB also forced my hand. When I save
% images with my saveAllFigs function, MATLAB turns all white annotations
% (text and lines) from white to black, rendering my scale bar useless.
invertedImage = imcomplement(cellImage);

% crop image according to the polygon's mirror calibration
croppedImage = imcrop(invertedImage, [1,100,1376,873]);
figure('name', strcat(fileName, '_', analysisDate, ' - psc_vs_light_polygon - Alexa image grid'));
hold on;
imshow(croppedImage, 'Border', 'tight');

% add polygon grid to figure
for row = 1:gridRows+1
    ypos = (row-1)*ymaxImg/gridRows;
    yline(ypos, 'Color', 'k', 'LineWidth', 0.5);
end
for col = 1:gridColumns+1
    xpos = (col-1)*xmaxImg/gridColumns;
    xline(xpos, 'Color', 'k', 'LineWidth', 0.5);
end

% add scale bar to figure
% ASSUMPTION ALERT: same parameters as DIC image
% line([x x], [y y])
line([scaleDistanceFromEdge scaleDistanceFromEdge+scaleBarSizeInPixels],...
    [ymaxImg-scaleDistanceFromEdge ymaxImg-scaleDistanceFromEdge],...
    'LineWidth', 2, 'Color', 'k');
% text(x, y, string)
text(scaleDistanceFromEdge,...
    ymaxImg-2*scaleDistanceFromEdge,...
    strcat(num2str(scaleBarSize), " μm"),...
    'FontSize', 10, 'Color', 'k');
hold off;

% re-size 
set(gcf,'InnerPosition',[2*innerWidth maxHeight-innerHeight innerWidth innerHeight]);


%% PLOT 3 - Polygon Heatmap 2 ================================================================================
% Made to be overlayed on top of cell image from rig

% re-organize data as a grid for heatmap
% the " .' " at the end makes sure that the sweeps are placed in the
% correct place in the grid, given an ordered polygon design - aka in each
% sweep the light moves one ROI to the right, or to the first ROI in
% the next row when it reached the end of the current row.
dataForHeatmap = reshape(sdFromPreLightHzPerROI,gridColumns,[]).';

% resize the heatmap matrix to match the size of the cropped image from the rig
% the original heatmap matrix will only be as big as the number of ROIs
% used in the polygon design. But the image from the rig will be a lot
% bigger (the size of the matrix is the size of image in pixels)
resizedHeatmap = imresize(dataForHeatmap, [size(croppedImage,1) size(croppedImage,2)], 'nearest');

% make heatmap without the heatmap function
figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - heatmap'));
imshow(resizedHeatmap,'Colormap',flipud(customColorMapPink),'DisplayRange', [heatmapMin,heatmapMax], 'Border', 'tight');
set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
c = colorbar;
c.Label.String = 'SDs from baseline firing rate';


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
for ROI = 1:totalROIs
    nexttile
    hold on;
    
    % convert cell array into regular matrix (for improved code
    % readability)
    relativeSweepsInThisROI = cell2mat(relativeSweepsInROI(ROI,:));
    
    for sweep = relativeSweepsInThisROI
        plot(cell2mat(tsBySweep(sweep)), ones(size(cell2mat(sweepNumberArrayBySweep(sweep))))+sweep, '|', 'Color', 'k')        
    end
    
    % zooming in and beautifying raster plot
    % range that works well for 3 sweeps
    axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur relativeSweepsInThisROI(1) relativeSweepsInThisROI(end)+2])
%     % range that works well for 10 sweeps
%     axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur relativeSweepsInThisROI(1)-20 relativeSweepsInThisROI(end)+20])

%     ylabel(strcat('Sweeps (', num2str(plannedSweepsPerROI), ')'));
%     xlabel('Time (s)');
    yticks([]);
    xticks(0:1:10);

    % adding light stim
    rectangle('Position', [lightOnsetTime relativeSweepsInThisROI(1)-20 lightDur relativeSweepsInThisROI(end)+100], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');

    % stop plotting things on this subplot
    hold off;
    
    % remove x labels from all plots except the last
    if ROI ~= gridRows * gridColumns
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
ylabel(t,strcat('Sweeps (', num2str(plannedSweepsPerROI), ')'));

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
for ROI = 1:totalROIs
    tsByROI = [];
    nexttile
    hold on;
    
    for sweep = cell2mat(relativeSweepsInROI(ROI,:))
        tsByROI = [tsByROI; cell2mat(tsBySweep(sweep))]; 
    end
    
    [Nroi, edges] = histcounts(tsByROI,edges);
    firingHzPerROI = Nroi/sweepsPerROI(ROI);
    histogram('BinEdges', 0:sweepSizeSec, 'BinCounts', firingHzPerROI, 'DisplayStyle', 'stairs', 'EdgeColor', 'k'); 

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
    if ROI ~= 1 
        yticklabels([]);
    end
    
    % remove x labels from all plots except the last
    if ROI ~= gridRows * gridColumns
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
    'ROI', ...
    'gridColumns', ...
    'gridRows', ...
    'isOrderedGrid', ...
    'plannedSweepsPerROI', ...    
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


% store key ROI-by-ROI data - each row is an ROI
filename = strcat(fileName, '_', analysisDate, " - firing_vs_light_polygon - ROI_by_ROI");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(dataCellMultipleRows);
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
    {'mouse', ...
    'date', ...
    'firstSweep', ...
    'lastSweep', ...
    'nSweeps', ...
    'gridColumns', ...
    'gridRows', ...
    'isOrderedGrid', ...
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
    'ROI', ...
    'sweepsPerROI', ...
    'preLight(Hz)AVG', ...
    'duringLight(Hz)AVG', ...
    'postLight(Hz)AVG', ...
    'preLight(Hz)STD', ...
    'duringLight(Hz)STD', ...
    'postLight(Hz)STD', ...
    'lightEffect(-1,0,1)', ...
    'medianSdFromPreLightHz'});
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the ROI_by_ROI xls file')
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


% end


%% Outdated code

% %% PLOT 2 - Polygon Heatmap 1 ================================================================================
% 
% % re-organize data as a grid for heatmap
% % the " .' " at the end makes sure that the sweeps are placed in the
% % correct place in the grid, given an ordered polygon design - aka in each
% % sweep the light moves one square to the right, or to the first square in
% % the next row when it reached the end of the current row.
% dataForHeatmap = reshape(dataPerROI(:,end),gridColumns,[]).';
% 
% % plot heatmap as an actual heatmap
% figure('name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - heatmap 1'));
% 
% % set the ColorLimits for consistent color-coding accross cells/mice
% % set the ColorMap to a flipped map with "flipud", so that more negative
% % values will be "hotter" and values closer to zero will be "colder".
% % remove labels from each square with CellLabelColor = 'none'.
% h = heatmap(dataForHeatmap, 'ColorLimits', [heatmapMin,heatmapMax], 'ColorMap', flipud(parula), 'CellLabelColor', 'none');
% 
% % removing labels from x axis
% cdl = h.XDisplayLabels;                                    % Current Display Labels
% h.XDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
% 
% % removing labels from y axis
% cdl = h.YDisplayLabels;                                    % Current Display Labels
% h.YDisplayLabels = repmat(' ',size(cdl,1), size(cdl,2));   % Blank Display Labels
% 
% % re-size 
% set(gcf,'InnerPosition',[2*innerWidth maxHeight-innerHeight innerWidth innerHeight]);
% 
% 
% %% PLOT 4 - Make a blended image with transparency ===========================================================
% % White where there was no inhibition of firing
% % Transparent where there was inhibition of firing beyond significance
% % criteria (firing decreased by -2 SDs from mean baseline firing rate)
% % Note that actual criteria is determined by user input "heatmapMin"
% % If heatmapMin = -2, then the criteria is "firing decreased by -2 SDs from
% % mean baseline firing rate"
% 
% % pilot code that is now obsolete:
% % figure('Name', '1')
% % background = im2double(croppedImage);
% % blendedImage = 0.5*resizedHeatmap + background;
% % blendedImage = im2uint8(blendedImage);
% % imshow(blendedImage)
% 
% % creating a transparency mask
% % trasform all negative values into positive ones: -resizedHeatmap
% transparencyMask = -resizedHeatmap;
% % default values are heatmapMax=0 and heatmapMin=-2
% % so default results are transparencyMaskMin=0 and transparencyMaskMax=2
% transparencyMaskMin = -heatmapMax;
% transparencyMaskMax = -heatmapMin;
% % anything between 0 and 2 --> normalized between 1 and 0 --> transparency gradient
% transparencyMask(transparencyMask>transparencyMaskMin & transparencyMask<transparencyMaskMax) = 1 - transparencyMask(transparencyMask>transparencyMaskMin & transparencyMask<transparencyMaskMax)/transparencyMaskMax;
% % anything below 0 --> 1 --> min transparency --> white
% transparencyMask(transparencyMask<=transparencyMaskMin) = 1;
% % anything above 2 --> 0 --> max transparency --> transparent
% transparencyMask(transparencyMask>=transparencyMaskMax) = 0;
% % sanity check that transparencyMask is between 0 and 1
% max(max(transparencyMask));
% min(min(transparencyMask));
% 
% figure('Name', strcat(fileName, '_', analysisDate, ' - firing_vs_light_polygon - transparent heatmap'))
% background = imshow(croppedImage, 'Border', 'tight');
% % make a white image the same size as the croppedImage
% whiteImage = cat(3, ones(size(croppedImage)), ones(size(croppedImage)), ones(size(croppedImage)));
% % fyi the placement of hold on and hold off in this code matter a lot for
% % whatever reason. Any other placement leads to image overwritting instead
% % of overlaying.
% hold on;
% % overlay white mask
% whiteMask = imshow(whiteImage);
% hold off;
% % apply transparency mask to white mask
% set(whiteMask, 'AlphaData', transparencyMask);
% 
% % set figure size to the same as the cropped cell image
% % FYI this only adjusts the outer figure size, not the inner figure size...
% set(gcf,'InnerPosition',[0 maxHeight-2.5*innerHeight innerWidth innerHeight]);
% % set(gcf,'OuterPosition',[1 1 outerWidth outerHeight]);



% % Keeping original variable names but assigning them to the new ones
% relativeSweepNumberPerSquare = relativeSweepsInROI;



% to properly plot histograms later, I need to keep track of sweeps per
% square taking into consideration any discarded sweeps. I will call this
% variable sweepsPerSquarePerSquare
% sweepsPerSquarePerSquareAll = sweepsPerROI;


% commented out legacy code
% % checking for incomplete sweeps and not analyzing incomplete sweeps - to
% % avoid this error: "Index exceeds the number of array elements (0)". 
% if numel(fieldnames(obj.sweeps)) <= obj.header.NSweepsPerRun  
%     lastSweepNumber = firstSweepNumber + numel(fieldnames(obj.sweeps)) -1 -discardedSweepsFromEnd;
%     allSweeps = firstSweepNumber:lastSweepNumber;
% end

% commented out legacy code
% % now checking if the total number of sweeps is a multiple of the total
% % number of squares. If not, remove the extra sweeps from the end.
% totalSweeps = length(allSweeps);
% totalSquares = gridColumns * gridRows;
% sweepsPerSquare = totalSweeps / totalSquares;
% if mod(totalSweeps, totalSquares) ~= 0
%     lastSweepNumber = lastSweepNumber - mod(totalSweeps, totalSquares);
%     allSweeps = firstSweepNumber:lastSweepNumber;
%     totalSweeps = length(allSweeps);
%     sweepsPerSquare = totalSweeps / totalSquares;
% end   


% % % store cell data with square-by-square info in a single row
% % filename = strcat(fileName, '_', analysisDate, " - firing_vs_light_polygon - square_by_square cell");
% % fulldirectory = strcat(savefileto,'\',filename,'.xls');        
% % dataInCellFormat = {};
% % dataInCellFormat = num2cell(dataCellSingleRow);
% % dataInCellFormat = [dataInCellFormat, cellImageFileNameDIC];
% % 
% % % create cell array with strings for naming square-by-square data
% % squareVariables = cell(1,size(dataPerROI,1)*size(dataPerROI,2));
% % for ROI = [1:totalROIs]
% %     startingIndex = ROI*9 - 8;
% %     squareVariables(startingIndex) = {strcat('sq', num2str(ROI), '-plannedSweepsPerROI')};
% %     squareVariables(startingIndex+1) = {strcat('sq', num2str(ROI), '-preLight(Hz)AVG')};
% %     squareVariables(startingIndex+2) = {strcat('sq', num2str(ROI), '-duringLight(Hz)AVG')};
% %     squareVariables(startingIndex+3) = {strcat('sq', num2str(ROI), '-postLight(Hz)AVG')};
% %     squareVariables(startingIndex+4) = {strcat('sq', num2str(ROI), '-preLight(Hz)STD')};
% %     squareVariables(startingIndex+5) = {strcat('sq', num2str(ROI), '-duringLight(Hz)STD')};
% %     squareVariables(startingIndex+6) = {strcat('sq', num2str(ROI), '-postLight(Hz)STD')};
% %     squareVariables(startingIndex+7) = {strcat('sq', num2str(ROI), '-lightEffect(-1,0,1)')};
% %     squareVariables(startingIndex+8) = {strcat('sq', num2str(ROI), '-medianSdFromPreLightHz')};
% % end
% % 
% % % label stored data
% % labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
% %     {'mouse', ...
% %     'date', ...
% %     'firstSweep', ...
% %     'lastSweep', ...
% %     'nSweeps', ...
% %     'gridColumns', ...
% %     'gridRows', ...
% %     'isOrderedGrid', ...
% %     'peaks(1)OrValleys(-1)', ...
% %     'highpassThreshold', ...
% %     'lowpassThreshold', ...
% %     'minPeakHeight', ...        
% %     'minPeakDistance', ...
% %     'lightExtensionFactor', ...
% %     'lightChannel', ...
% %     'singleLightPulse(1)orTrain(0)', ...
% %     'heatmapMin', ...
% %     'heatmapMax', ...
% %     'lightPulseDur(s)', ... 
% %     'lightStimFreq(Hz)', ...
% %     'lightDur(s)', ...   
% %     'APhalfWidth(ms)', ...
% %     'APduration(ms)', ...
% %     'baselineHzMean', ...
% %     'baselineHzSD', ...
% %     'baselineHzCV', ...
% %     'baselineISImean', ...
% %     'baselineISIsd', ...
% %     'baselineISIcv', ...
% %     'isDA(0,1)', ...
% %     'isIrregular(0,1)', ...
% %     squareVariables{:}, ...
% %     'cellImageFileNameDIC'});
% % writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
% % disp('I saved the square_by_square cell xls file')
