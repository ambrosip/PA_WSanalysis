%{ 
DOCUMENTATION
Created: 2022 09 21
Edited on: 2022 11 04; 2022 12 24; 2022 12 25; 2024 10 29
Major github merge: 2022 12 24
Last edits:
    2024 10 - removed lines about images
    2024 11 - turned some user inputs into function variables for batch
    processing of data 
Works? Yes
Author: PA

This function is used to analyze light-evoked PSCs in sCRACM experiments
using the polygon. This function was derived from psc_vs_light_polygon, but
went through A LOT of rewritting.
The data organization is very different, to better handle discarded sweeps
without losing track of which sweeps belong to which ROIs (Regions Of
Interest - subarea of the field that was illuminated).
The code was further improved to accomodate custom ordered grids in
addition to ordered grids.
The code was also modified to handle "problem files" in which the timing of
o-stims was messed up during acquisition.
FYI I sometimes use the words "square" and "ROI" interchangeably.
2022 12 25 edit: custom crop to accomodate changed in polygon's mirror
calibration
By 2022 12 25, these are the channel standards:
ch1: data
ch2: voltage command
ch3: LED command
ch4: polygon command

This function uses a WSfile object. Any changes to the name of this
function and/or its inputs (the stuff within parenthesis) need to be added
to method section of "WSfile.m"

Adaptation tree:
    psc_vs_light_polygon

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
        For problem files, these are the criteria for choosing discarded
        sweeps: (1) 2 voltage test pulses or 2 light pulses in the same
        sweeps; or (2) test/light pulses too close to the beginning or the
        end of the sweep; or (3) incorrect ROI/sweep assignment --> if you
        discarded a sweep due to multiple light pulses, the next sweeps are
        going to be incorrectly assigned to their ROIs.

    - lightChannel: channel where the info about the light stim is stored.
    Usually 2 (LED 1), 3 (LED 2), or 4 (polygon). If a digital output was
    used (my default for the polygon), the only usable info in this channel
    is timing (the start/end of the light pulse and duration of the light
    pulse). If an analog output was used (my default for the Mightex's
    LEDs), this channel has timing and amplitude information (0-100% LED
    power is represented as 0-5 V). When using small LED power, the
    signal/noise ratio for timing info is shitty, so I genrally opted to
    use the polygon channel for timing info and the LED channel for
    amplitude info (see variable below).

    - ledPowerChannel: channel where the info about the LED power is
    stored if the LED driver was triggered with an analog output (0-5 V).

    - singleLightPulse: boolean variable (0 or 1) determining whether this
    sweep has a single light pulse (1) or a train of light pulses (0).
    Added to avoid errors

    - inwardORoutward: tells MATLAB to look for an inward current (-1) or
    an outward current (+1)

    - baselineDurationInSeconds: tells MATLAB where to calculate the
    baseline for baseline subtraction

    - lightPulseAnalysisWindowInSeconds: how long after light onset to look
    for light evoked responses

    - thresholdInDataPts: min number of points changing in a monotonic
    fashion during the lightPulseAnalysisWindow to find PSC onset. Early in
    Ambrosi et al 2022, this was set to 10. Later, it was set to 5 (to
    acommodate the identification of smaller PSCs). In the polygon
    analysis, it was initially set to 5. However, I noticed that my code
    was overestimating the failure rate (oIPSCs were clearly detected by
    visual inspection, but failed to be detected by the code). This was due
    to noise in the recording hindering the detection of monotonic
    changes. To address this, I added an extra step in the analysis to
    smooth the data a little bit. I am using a moving average filter with a
    very narrow span (see smoothSpan) so that minor dataPoint to dataPoint
    fluctuations are dampened without significantly altering the kinetics
    of the recorded PSCs. After this change, setting thresholdInDataPts to 5
    is no longer good - it's too easy: the code underestimates the failure
    rate and find oIPSCs were there are none (based on visual inspection).
    Thus, the current default of this value is 10 --> a monotonic change in
    current for at least 1 ms.

    - amplitudeThreshold: min amplitude of light-evoked response to
    consider it a real response and not a failure. Do not need to worry
    about sign (+ or -). The default 25 was chosen based on the noise level of
    my rig, which is ~20 pA (as of 2022-09-21).

    - smoothSpan: number of data points averaged in the smoothing filter.
    If you change this value, you might have to change thresholdInDataPts
    as well. The default value is 3, which basically means that each
    dataPoint will be replaced by the average of 3 dataPoints. This mild
    smoothing ensures that noise fluctuations are dampened without lowpass
    filterding the data too much. If you increase smoothSpan, the kinetics
    of the detected PSCs will be slower than the real kinetics. A quick
    test revealed that setting smoothSpan to 5 is already enough to mess
    with onset latency and peak amplitude a bit.

    - discardROIsWithLowFreq: boolean variable (0 or 1) determining whether
    to discard "PerROI" data from ROIs with success rate < 50% (1) or not
    (0).

    - problemFile: boolean variable (0 or 1) determining whether the timing
    of the o-stims got messed up during acquisition, which will require the
    code to look for the unique timing of each o-stim in all sweeps and
    discard sweeps in which a huge fucking mess would occur.
        For problem files, these are the criteria for choosing discarded
        sweeps: (1) 2 voltage test pulses or 2 light pulses in the same
        sweeps; or (2) test/light pulses too close to the beginning or the
        end of the sweep; or (3) incorrect ROI/sweep assignment --> if you
        discarded a sweep due to multiple light pulses, the next sweeps are
        going to be incorrectly assigned to their ROIs.

    - rsTestPulseOnsetTime: onset of voltage step used to calculate series
    resistance.

    - autoRsOnsetTime: boolean variable (0 or 1) determining whether to use
    rsTestPulseOnsetTime to calculate Rs (0) or to look for exact location
    of the voltage pulse based on the voltageCmdChannel (1).

    - voltageCmdChannel: channel where the info about the voltage command
    is stored. Usually 2.

    - ymin: for illustration purposes, use this value as the min range for
    the y-axis when plotting current vs time.

    - ymax: for illustration purposes, use this value as the max range for
    the y-axis when plotting current vs time.

    - cellImageFileNameDIC: file containing DIC image for analyzed cell.
 
    - cellImageFileNameAlexa: file containing Alexa image for analyzed cell

    - cellImageDir: folder where cell image file is located

    - leftCrop: the polygon mirror area is smaller than the field of view
    of the rig camera. To properly overlay the grid onto the images taken
    from the rig, I need to crop the image from ocular. This crop requires
    four parameters: leftCrop, rightCrop, topCrop and bottomCrop. The first
    corresponds to the number of pixels in the x-axis betwween the left-most 
    edge of the photo and the left-most edge of the mirror grid.

    - rightCrop: number of pixels between the right-most edge of the
    polygon's mirror grid and the right edge of the ocular picture.

    - topCrop: number of pixels between the top-most edge of the
    polygon's mirror grid and the top edge of the ocular picture.

    - bottomCrop: number of pixels between the bottom edge of the
    polygon's mirror grid and the bottom edge of the ocular picture.

    - savefileto: save CSV files to this folder.

INPUTS defaults:
    % Affects data analysis - Organizing data by o-stim grid
    gridColumns = 5;
    gridRows = 5;
    orderedGrid = 0;    % 0 if NOT ordered, 1 if ordered  
    orderOfROIs = [8 16 14 23 3 10 12 25 21 19 6 18 1 5 11 4 22 15 13 7 2 20 24 17 9]';     % this is the order of the design 5x5 spaced out

    % Affects data analysis - Finding/quantifyting oIPSCs
    discardedSweeps = [];
    discardedSweepsFromEnd = 0;
    lightChannel = 4;
    ledPowerChannel = 3;
    singleLightPulse = 1; 
    inwardORoutward = -1;               
    baselineDurationInSeconds = 0.01;
    lightPulseAnalysisWindowInSeconds = 0.02;
    thresholdInDataPts = 5;     
    amplitudeThreshold = 25;
    smoothSpan = 3; 
    discardROIsWithLowFreq = 1; 
    problemFile = 0;

    % Affects data analysis - Calculating Rs
    rsTestPulseOnsetTime = 1;
    autoRsOnsetTime = 0;
    voltageCmdChannel = 2;

    % Affects data display - oIPSC amplitude range (in pA): 
    ymin = -1500;           %-2050      -3600      -1500
    ymax = 600;             %50         600        600

    % Affects data display - polygon grid overlay & crop
    cellImageFileNameDIC = 's3c1_dic_1x.tif';
    cellImageFileNameAlexa = 's3c1_MAX_Stack Rendered Paths.tif';
    cellImageDir = 'D:\NU server\Priscilla - BACKUP 20200319\Ephys\2022\20220726 m000 dls loop';
    leftCrop = 0;     % in pixels              % ALERT! this is new (2022-12-25)
    rightCrop = 0;    % in pixels
    topCrop = 100;    % in pixels
    bottomCrop = 51;  % in pixels

    % Affects data saving:
    savefileto = 'D:\Temp\From MATLAB 2022 08 31 psc';    
    
OUTPUTS:
    Fig
    XLS

ASSUMPTIONS: 
    - The same ROI order was used throughout the file
    - Width and frequency of polygon stim is the the same as LED o-stim
    (only relevant if you use the polygon channel as the light stim
    channel)

BEWARE:
    - ctrl+f for ALERT and ASSUMPTION

TO DO:
    - noticed weird bug on 8/31/22: if I include the last sweep on the
    discarded sweeps list, matlab does not plot any of the tiled niceplots!
    - check github merge consequences
    - add plot of avg traces (accross sweeps) per ROI & avg oIPSC peak per
    ROI ignoring failures (similar analysis pipeline from Mao Lab and
    Gordon Sheppard)
%}

% obj = m729.s0246;
% obj = m729.s0111;
% obj = m729.s0371;
% obj = m731.s0007

function psc_vs_light_polygon_batch(obj, leftCrop, rightCrop, topCrop, bottomCrop, gridColumns, gridRows, orderOfROIs, discardedSweeps, saveDir, addNamePrefix, namePrefix)
%%  USER INPUT ============================================================

orderedGrid = 0;       % 0 if NOT ordered, 1 if ordered

% Affects data analysis - Finding/quantifyting oIPSCs
lightChannel = 4;
ledPowerChannel = 3;
singleLightPulse = 1; 
inwardORoutward = -1;                       % 1 (positive) is outward; -1 (negative) in inward
baselineDurationInSeconds = 0.01;
lightPulseAnalysisWindowInSeconds = 0.015;  % ALERT: changed from 0.02 to 0.01 to 0.015 on 2022-09-21. Changed to 0.045 on 2024-10-29 - matches the hard-coded xmax used in the tiled niceplots
thresholdInDataPts = 8;                     % ALERT! Changed from 10 to 5 to 10 to 8 (2022-09-24)
amplitudeThreshold = 25;                    % ALERT! this is new (2022-09-21)
smoothSpan = 3;                             % ALERT! this is new (2022-09-23)
discardROIsWithLowFreq = 0;                 % ALERT! this is new (2022-09-24). Changed from 1 to 0 on 2025-02-05
problemFile = 0;                            % ALERT! this is new (2022-09-24)

% Affects data analysis - Calculating Rs
rsTestPulseOnsetTime = 1;
autoRsOnsetTime = 1;
voltageCmdChannel = 2;

% Affects data display - oIPSC amplitude range (in pA):
% use the [-2050 50] range for a scale bar of 300 pA
% use the [-3600 600] range for a scale bar of 600 pA
ymin = -2050;           %-2050      -3600       
ymax = 50;             %50         600        

% Affects data saving:
% savefileto = 'Z:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries\2023\2023 09 18 scracm asc m928';
savefileto = saveDir;

% Affects data display
gridFillHorizontal = 1;     % 0.067 for 15 column
gridFillVertical = 1;       % 0.067 for 15 line


%% PREP - get monitor info for plot display organization =====================================================

monitorPositions = get(0, 'MonitorPositions' );
maxHeight = monitorPositions(1,4) - 100;
maxWidth = monitorPositions(1,3) - 100;


%% PREP - get info about cropped image size

cropXmin = 1 + leftCrop;
cropWidth = 1376 - rightCrop - leftCrop;
cropYmin = 1 + topCrop;
cropHeight = 1024 - topCrop - bottomCrop;


%% PREP - get info from h5 file and create arrays ============================================================

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% getting info from h5 file
mouseNumber = getMouseNumber(obj);
experimentDate = getExperimentDate(obj);
samplingFrequency = obj.header.Acquisition.SampleRate;
sweepDurationInSec = obj.header.SweepDuration;
sweepDurationInDataPts = sweepDurationInSec * samplingFrequency;

% get ingo from h5 file or use user-defined name prefix
if addNamePrefix == 1
    fileName = namePrefix;
else
    fileName = obj.file;
end

% calculating variables based on user input - the Rs-related variables
% might change later in the code depending on the value in autoRsOnsetTime
baselineDurationInDataPts = baselineDurationInSeconds * samplingFrequency;
lightPulseAnalysisWindowInDataPts = lightPulseAnalysisWindowInSeconds * samplingFrequency;
rsBaselineDataPointInterval = ((rsTestPulseOnsetTime-0.05)*samplingFrequency):(rsTestPulseOnsetTime*samplingFrequency);
rsFirstTransientDataPointInterval = (rsTestPulseOnsetTime*samplingFrequency):(rsTestPulseOnsetTime+0.0025)*samplingFrequency;

% getting sweep numbers from file name
[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj); 

% ALERT: NEED TO TEST
% checking for incomplete sweeps and adding them to the list of discarded sweeps
% to avoid this error: "Index exceeds the number of array elements (0)". 
if numel(fieldnames(obj.sweeps)) < obj.header.NSweepsPerRun  
    lastCompleteSweep = firstSweepNumber + numel(fieldnames(obj.sweeps)) - 1;
    discardedSweeps = [discardedSweeps, lastCompleteSweep:lastSweepNumber];
end

% ALERT: NEED TO TEST
% checking for round number of sweeps per ROI
% to avoid this error: "Size arguments must be real integers".
totalROIs = gridColumns * gridRows;
extraSweeps = rem(length(allSweeps),totalROIs);
if extraSweeps ~= 0
    lastCompleteSweep = lastSweepNumber - extraSweeps;
    discardedSweeps = [discardedSweeps, lastCompleteSweep+1:lastSweepNumber];
end

% add 'subset' to fileName in case discardedSweeps is not empty
% to prevent overwritting of saved files with the whole dataset
if ~isempty(discardedSweeps)
    fileName = strcat(fileName, 'subset_');
    allSweeps = setdiff(allSweeps, discardedSweeps);
end

% organizing sweeps according to total ROIs and discarding discardedSweeps
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

% creating matrixes/arrays that will be filled later
allRs = [];
yBaselineSubAll = [];
allChargesPerSweep = [];
allLightEvokedCurrentsPerSweep = [];
allFailureAssessmentsPerSweep = [];
allLightEvokedResponseOnsetLatencyInMilliSecondsPerSweep = [];
allLightEvokedResponsePeakLatencyInMilliSecondsPerSweep = [];
allTimeTo10percentOfPeakInMilliSecondsPerSweep = [];
allTimeTo90percentOfPeakInMilliSecondsPerSweep = [];
allRiseTimeInMilliSecondsPerSweep = [];
baselineCurrentAll = [];
data = [];
dataInROI = cell(totalROIs,1);
yAroundLightPulseAll = [];
heavyEditingAll = [];

% clearing cell arrays for safety
% (if you use this code as a script, data will not linger from one analysis
% to the next)
meanDataInROI = [];
chargeInROI = {};
peaksInROI = {};
onsetLatenciesInROI = {};
riseTimesInROI = {};
failuresInROI = {};
meanChargePerROI = [];
meanPeakPerROI = [];
meanOnsetLatencyPerROI = [];
meanRiseTimePerROI = [];
failureRatePerROI = [];
stdOnsetLatencyPerROI = [];


%% ROI BY ROI AND SWEEP BY SWEEP ANALYSIS ===================================================================================

% get data from sweeps in file (only the subset we will analyze)
% ROI: region on interest (subarea illuminated)
for ROI = 1:totalROIs  
    
    column = 0;
    
    for sweepNumber = cell2mat(sweepsInROI(ROI))
        
        column = column + 1;
           
        % get light stim data
        [xLight,yLight] = obj.xy(sweepNumber, lightChannel); 

        % get light stim parameters
        % if you are using a TTL pulse (5V) to control the LED, use code #1
        % if you are using an analog output (0-5V), use code #2
        % ALERT let's just use the polygon TTL pulse for now, since the analog
        % output is soooo small, making the detection of the light pulse start
        % and end times really hard with simple methods.
        % aka, use code #1 and change the lightChannel to the polygon channel
        % instead of the LED channel.

        % code #1 - works
        % look for a big change 
        % I was staring at this code and it does not make sense
        % yLight>1 is boolean bruh - it tells you whether yLight is > 1 for each data point - that is useless
%         lightPulseStart = find(diff(yLight>1)>0);
%         lightPulseEnd = find(diff(yLight<1)>0);

        % code #1 that makes more sense
        % look for a very positive and very negative derivative
        % ALERT: have not tested this code for multiple light pulses
        lightPulseStart = find(diff(yLight)>1);
        lightPulseEnd = find(diff(yLight)<-1);

%         % code #2 - does not always work
%         % look for a small change
%         % to avoid artifacts, use both the derivative and the absolute value of
%         % ych2 to find the start and end times of each light pulse
%         % ALERT need to check this code with a train o-stim
%         lightPulseStart = find(diff(ych2)>0.075 & ych2(1:end-1)<0.05);
%         lightPulseEnd = find(diff(ych2)<-0.075 & ych2(1:end-1)>0.05);

        % continue to get light stim info
        lightOnsetTime = lightPulseStart(1)/samplingFrequency;                       % in seconds
        stimDur = (lightPulseEnd(end)-lightPulseStart(end))/samplingFrequency;       % duration of each light pulse in the train (s)

        % if the light stim is a train (singleLightPulse = 0), compute light
        % train information. 
        % also set xmax for plotting later (TO DO: test this xmax)
        if singleLightPulse == 0
            stimInterval = (lightPulseStart(2)-lightPulseStart(1))/samplingFrequency;    % interval between each pulse (s)
            stimFreq = 1/stimInterval;                                                   % frequency of the light stim (Hz)
            lightDur = (lightPulseStart(end)-lightPulseStart(1))/samplingFrequency + stimInterval;  % duration of the whole light train stim (s) 
            xmax = lightOnsetTime + lightDur + baselineDurationInSeconds;

        % if the light stim is a single pulse (singleLightPulse = 1), set the
        % train information to the following values (to avoid errors)
        % also set xmax for plotting later
        else
            stimInterval = 0;
            stimFreq = 1;
            lightDur = stimDur;
            % ALERT: hard-coded xmax
            % jesus christ I had to change this xmax from data-based to
            % kind of hard-coded so that the time scale bar has a round number.
            % might be worth messing with the scale bar in the future.
            % commented out below are the previous iterations of xmax
            % calculations:
            xmax = lightOnsetTime + lightDur + 0.04;
%             xmax = lightOnsetTime + lightDur + 2*baselineDurationInSeconds; 
%             xmax = lightOnsetTime + lightDur + baselineDurationInSeconds; 
%             xmax = lightOnsetTime+0.2;        
        end    

        % set xmin for plotting later
        xmin = lightOnsetTime - baselineDurationInSeconds;
        %---------------------------------------------------------------------- 
        
        % get LED power
        % since we're using the polygon channel as the light stim channel,
        % I'm just gonna have to get data from the actual LED channel for
        % power measurement 
        [xLED,yLED] = obj.xy(sweepNumber, ledPowerChannel);
        ledPower = round(max(yLED),1);
        %----------------------------------------------------------------------
        
        % get data from channel 1 (current recording)
        [x,y] = obj.xy(sweepNumber, 1);
              
        % checking for problem sweeps in which the total duration of the
        % recorded data does not match the planned duration
        heavyEditing = 0;
        plannedSweepDurationInDataPoints = obj.header.SweepDuration * samplingFrequency;
        if size(x,1) > plannedSweepDurationInDataPoints
            % remove extra data to avoid errors in the concatenation
            % happening next
            heavyEditing = 1;
            heavyEditingAll = [heavyEditingAll, heavyEditing];
            x(plannedSweepDurationInDataPoints+1:end) = [];
            y(plannedSweepDurationInDataPoints+1:end) = [];
        end
        %----------------------------------------------------------------------      
        
        % smooth data with moving average filter
        y = smooth(y,smoothSpan);

        % baseline subtraction
        baselineStart = lightPulseStart(1) - baselineDurationInDataPts;
        yBaselineSub = y-mean(y(baselineStart:lightPulseStart(1)));
        baselineCurrent = mean(y(baselineStart:lightPulseStart(1)));
        %----------------------------------------------------------------------       
        
        % if the timing of the light pulse changes accross sweeps, let's
        % keep only the data around the lightpulse - so that we plot the
        % right thing later
        % saving a subset of the data around the lightpulse
        xminDataPts = round(xmin * samplingFrequency);
        xmaxDataPts = round(xmax * samplingFrequency);
        xAroundLightPulse = x(xminDataPts:xmaxDataPts);
        yAroundLightPulse = yBaselineSub(xminDataPts:xmaxDataPts);
        yAroundLightPulseAll = [yAroundLightPulseAll, yAroundLightPulse]; 
        %----------------------------------------------------------------------             
               
        % saving data for niceplot
        % y data for each sweep is in a column
        yBaselineSubAll = [yBaselineSubAll, yBaselineSub]; 
        baselineCurrentAll = [baselineCurrentAll, baselineCurrent];
        %----------------------------------------------------------------------    

        % checking relevant timepoints for series resistance calculation
        % if user indicated auto Rs check, that means that the timing of the Rs
        % test pulse changes from sweep to sweep. So we need to calculate
        % rsTestPulseOnsetTime based on the data on voltageCmdChannel
        if autoRsOnsetTime == 1
            [xV,yV] = obj.xy(sweepNumber, voltageCmdChannel);
            rsTestPulseDataPoint = find(diff(yV)<-1);
            rsTestPulseOnsetTime = rsTestPulseDataPoint/samplingFrequency;
            rsBaselineDataPointInterval = ((rsTestPulseOnsetTime-0.05)*samplingFrequency):(rsTestPulseOnsetTime*samplingFrequency);
            rsFirstTransientDataPointInterval = (rsTestPulseOnsetTime*samplingFrequency):(rsTestPulseOnsetTime+0.0025)*samplingFrequency;
            
            % trying to troubleshoot an error that I don't understand
            rsBaselineDataPointInterval = round(rsBaselineDataPointInterval(1)):round(rsBaselineDataPointInterval(end));
            rsFirstTransientDataPointInterval = round(rsFirstTransientDataPointInterval(1)):round(rsFirstTransientDataPointInterval(end));
        end
        %----------------------------------------------------------------------
        
        % calculating series resistance
        rsBaselineCurrent = mean(y(rsBaselineDataPointInterval));  
        rsTransientCurrent = min(y(rsFirstTransientDataPointInterval));
        dCurrent = rsTransientCurrent-rsBaselineCurrent;
        dVoltage = -5;  % ASSUMPTION ALERT
        seriesResistance = 1000*dVoltage/dCurrent; %mV/pA equals Gohm

        % put series resistance values from each sweep into a different column
        allRs = [allRs, seriesResistance];
        %----------------------------------------------------------------------
        
        % here comes the actual psc_vs_light analysis:
        % clean up matrices that will be used in the next loop
        % "all" in the beginning refers to all light pulses in a sweep
        allChargesPerLightPulse = [];
        allLightEvokedCurrentsAmp = [];
        allLightEvokedCurrentsLoc = [];
        allFailureAssessments = [];
        allLightEvokedResponseOnsetLatencyInMilliSeconds = [];
        allLightEvokedResponsePeakLatencyInMilliSeconds = [];
        allTimeTo10percentOfPeakInMilliSeconds = [];
        allTimeTo90percentOfPeakInMilliSeconds = [];
        allRiseTimeInMilliSeconds = [];

        % loop through all light pulses in the train
        % this code is overkill if there is only one light pulse but whatever
        for pulseOnset = lightPulseStart.'
            
            % this is used to find light-evoked current onset
            afterLightDataPoint = pulseOnset + lightPulseAnalysisWindowInDataPts; 

            % get charge
            dataSubsetForCharge = yBaselineSub(pulseOnset:afterLightDataPoint);
            chargeAfterLightPulse = trapz(dataSubsetForCharge,1); % pA x data point
            chargeAfterLightPulse = chargeAfterLightPulse * 1/samplingFrequency; % pA x second

            % get onset latency of lightEvokedCurrents
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
            lightEvokedResponseOnsetLatencyInDataPoints = find(conv2(sign(diff(y(pulseOnset:afterLightDataPoint))), zeros(thresholdInDataPts,1)+(inwardORoutward), 'valid')==thresholdInDataPts, 1 );
           
            % if onset was detected, get amplitude and location (index) of lightEvokedCurrents
            % also calculate rise time
            if ~isempty(lightEvokedResponseOnsetLatencyInDataPoints)      
                % look for peak current AFTER onset to avoid misleading riseTimes
                onsetDataPoint = pulseOnset + lightEvokedResponseOnsetLatencyInDataPoints;
                afterOnsetDataPoint = onsetDataPoint + lightPulseAnalysisWindowInDataPts;                
                % outward current --> inwardORoutward is +1
                % inward current --> inwardORoutward is 0 or -1
                if inwardORoutward == 1 
                    [lightEvokedCurrentAmp, lightEvokedCurrentLoc] = max(yBaselineSub(onsetDataPoint:afterOnsetDataPoint));
                    timeTo10percentOfPeakInDataPoints = find(yBaselineSub(onsetDataPoint:afterOnsetDataPoint) >= 0.1 * lightEvokedCurrentAmp, 1);
                    timeTo90percentOfPeakInDataPoints = find(yBaselineSub(onsetDataPoint:afterOnsetDataPoint) >= 0.9 * lightEvokedCurrentAmp, 1);
                else
                    [lightEvokedCurrentAmp, lightEvokedCurrentLoc] = min(yBaselineSub(onsetDataPoint:afterOnsetDataPoint));
                    timeTo10percentOfPeakInDataPoints = find(yBaselineSub(onsetDataPoint:afterOnsetDataPoint) <= 0.1 * lightEvokedCurrentAmp, 1);
                    timeTo90percentOfPeakInDataPoints = find(yBaselineSub(onsetDataPoint:afterOnsetDataPoint) <= 0.9 * lightEvokedCurrentAmp, 1);
                end         
            % if onset was not detected, NaN it all! 
            % aka there is no light-evoked respone in response to this light pulse
            else
                lightEvokedCurrentAmp = NaN;
                lightEvokedCurrentLoc = NaN;
                timeTo10percentOfPeakInDataPoints = NaN;
                timeTo90percentOfPeakInDataPoints = NaN;
            end
            
            % use user-defined amplitude threshold to test if there was a
            % light-evoked response in this sweep or if we should consider
            % this a failure
            % ALERT: there are two criteria for defining failures:
            % 1) absence of fast current deflection within 15 ms of light
            % onset (implemented through onset latency calculation)
            % 2) peak current within 15 ms of deflection onset <
            % amplitudeThreshold (implemented through the code below)
            if abs(lightEvokedCurrentAmp) > amplitudeThreshold
                isFailure = 0;
            else
                isFailure = 1;
            end
            
            % Convert data points to milliseconds (1000 multiplication is conversion from seconds to milliseconds)
            lightEvokedResponseOnsetLatencyInMilliSeconds = 1000*lightEvokedResponseOnsetLatencyInDataPoints/samplingFrequency;
            lightEvokedResponsePeakLatencyInMilliSeconds = 1000*lightEvokedCurrentLoc/samplingFrequency;
            timeTo10percentOfPeakInMilliSeconds = 1000*timeTo10percentOfPeakInDataPoints/samplingFrequency;
            timeTo90percentOfPeakInMilliSeconds = 1000*timeTo90percentOfPeakInDataPoints/samplingFrequency;
            riseTimeInMilliSeconds = timeTo90percentOfPeakInMilliSeconds - timeTo10percentOfPeakInMilliSeconds;       

            % check these things to avoid errors while exporting data
            if isempty(lightEvokedResponseOnsetLatencyInMilliSeconds)
                lightEvokedResponseOnsetLatencyInMilliSeconds = NaN;
            end

            if isempty(lightEvokedResponsePeakLatencyInMilliSeconds)
                lightEvokedResponsePeakLatencyInMilliSeconds = NaN;
            end

            if isempty(timeTo10percentOfPeakInMilliSeconds)
                timeTo10percentOfPeakInMilliSeconds = NaN;
            end

            if isempty(timeTo90percentOfPeakInMilliSeconds)
                timeTo90percentOfPeakInMilliSeconds = NaN;
            end

            if isempty(riseTimeInMilliSeconds)
                riseTimeInMilliSeconds = NaN;
            end        

            % put all info from a particular sweep in a row (each column is a light pulse)
            allChargesPerLightPulse = [allChargesPerLightPulse, chargeAfterLightPulse];
            allLightEvokedCurrentsAmp = [allLightEvokedCurrentsAmp, lightEvokedCurrentAmp];
            allLightEvokedCurrentsLoc = [allLightEvokedCurrentsLoc, lightEvokedCurrentLoc];
            allFailureAssessments = [allFailureAssessments, isFailure];
            allLightEvokedResponseOnsetLatencyInMilliSeconds = [allLightEvokedResponseOnsetLatencyInMilliSeconds, lightEvokedResponseOnsetLatencyInMilliSeconds];        
            allLightEvokedResponsePeakLatencyInMilliSeconds = [allLightEvokedResponsePeakLatencyInMilliSeconds, lightEvokedResponsePeakLatencyInMilliSeconds];        
            allTimeTo10percentOfPeakInMilliSeconds = [allTimeTo10percentOfPeakInMilliSeconds, timeTo10percentOfPeakInMilliSeconds];
            allTimeTo90percentOfPeakInMilliSeconds = [allTimeTo90percentOfPeakInMilliSeconds, timeTo90percentOfPeakInMilliSeconds];
            allRiseTimeInMilliSeconds = [allRiseTimeInMilliSeconds, riseTimeInMilliSeconds];

        end

        % store info from all sweeps (each row is a sweep, each column is a light pulse)
        allChargesPerSweep = [allChargesPerSweep; allChargesPerLightPulse];
        allLightEvokedCurrentsPerSweep = [allLightEvokedCurrentsPerSweep; allLightEvokedCurrentsAmp];
        allFailureAssessmentsPerSweep = [allFailureAssessmentsPerSweep; allFailureAssessments];
        allLightEvokedResponseOnsetLatencyInMilliSecondsPerSweep = [allLightEvokedResponseOnsetLatencyInMilliSecondsPerSweep; allLightEvokedResponseOnsetLatencyInMilliSeconds];
        allLightEvokedResponsePeakLatencyInMilliSecondsPerSweep = [allLightEvokedResponsePeakLatencyInMilliSecondsPerSweep; allLightEvokedResponsePeakLatencyInMilliSeconds];
        allTimeTo10percentOfPeakInMilliSecondsPerSweep = [allTimeTo10percentOfPeakInMilliSecondsPerSweep; allTimeTo10percentOfPeakInMilliSeconds];
        allTimeTo90percentOfPeakInMilliSecondsPerSweep = [allTimeTo90percentOfPeakInMilliSecondsPerSweep; allTimeTo90percentOfPeakInMilliSeconds];
        allRiseTimeInMilliSecondsPerSweep = [allRiseTimeInMilliSecondsPerSweep; allRiseTimeInMilliSeconds];    
                   
        % store "raw" data in cell array
        % each ROI is stored in a row
        % each sweep within a ROI is stored in a column
        % adding an if statement to deal with problematic files
        % if any sweep went through heavyEditing due to inconsistent o-stim
        % timing, do this:
%         if any(heavyEditingAll)   % UGH this is not good enough cuz
%         problem sweeps might happen AFTER normal sweeps, messing up with
%         the data organization and giving errors.
        % let's go with a user-input based check UGH
        if problemFile
            dataInROI(ROI,column) = {yAroundLightPulse};
        else
            dataInROI(ROI,column) = {yBaselineSub};
        end
        
        % store key data for statistics (mean, SD, median) later
        chargeInROI(ROI,column) = {allChargesPerSweep};
        peaksInROI(ROI,column) = {allLightEvokedCurrentsAmp};
        onsetLatenciesInROI(ROI,column) = {allLightEvokedResponseOnsetLatencyInMilliSeconds};
        riseTimesInROI(ROI,column) = {allRiseTimeInMilliSeconds};
        failuresInROI(ROI,column) = {isFailure};
        
        % store key data that will be exported (each row is a sweep)
        data = [data; ...
            mouseNumber, ...
            experimentDate, ...
            sweepNumber, ...   
            ROI, ...
            gridColumns, ...
            gridRows, ...
            orderedGrid, ...
            plannedSweepsPerROI, ...
            lightChannel, ...
            ledPowerChannel, ...
            singleLightPulse, ...        
            inwardORoutward, ...
            baselineDurationInSeconds, ...
            lightPulseAnalysisWindowInSeconds, ...
            thresholdInDataPts, ...
            amplitudeThreshold, ...
            smoothSpan, ...
            problemFile, ...
            rsTestPulseOnsetTime, ...
            autoRsOnsetTime, ...
            voltageCmdChannel, ...
            stimDur, ... 
            stimFreq, ...
            lightDur, ...
            ledPower, ...
            seriesResistance, ...
            baselineCurrent, ...
            heavyEditing];         
    end 
    
    meanDataInROI = [meanDataInROI, mean(cell2mat(dataInROI(ROI,:)),2)];
    
    % ALERT: the code below might break down if you have multiple light
    % pulses, so I'm throwing a message out there
    if singleLightPulse ~= 1
        disp('bruh the "perROI" stuff was not made to handle multiple light pulses')
    end

    % meanChargePerROI(ROI) = mean(cell2mat(chargeInROI(ROI,:)), 'omitnan');
    meanPeakPerROI(ROI) = mean(cell2mat(peaksInROI(ROI,:)), 'omitnan');
    meanOnsetLatencyPerROI(ROI) = mean(cell2mat(onsetLatenciesInROI(ROI,:)), 'omitnan');    
    meanRiseTimePerROI(ROI) = mean(cell2mat(riseTimesInROI(ROI,:)), 'omitnan');   
    
    % nnz: number of non-zero elements (aka number of failures)
    failureRatePerROI(ROI) = 100*nnz(cell2mat(failuresInROI(ROI,:)))/sweepsPerROI(ROI);
       
    % if there are more than 2 sweeps with onset latency in this ROI, calculate SD
    if nnz(~isnan(cell2mat(onsetLatenciesInROI(ROI,:)))) > 2
        stdOnsetLatencyPerROI(ROI) = std(cell2mat(onsetLatenciesInROI(ROI,:)), 'omitnan');
    else
        stdOnsetLatencyPerROI(ROI) = NaN;
    end
            
end

% add pulse-by-pulse info to sweep-by-sweep data
data = [data, ... 
    allChargesPerSweep, ...
    allFailureAssessmentsPerSweep, ...
    allLightEvokedCurrentsPerSweep, ...
    allLightEvokedResponseOnsetLatencyInMilliSecondsPerSweep, ...
    allLightEvokedResponsePeakLatencyInMilliSecondsPerSweep, ...
    allTimeTo10percentOfPeakInMilliSecondsPerSweep, ...
    allTimeTo90percentOfPeakInMilliSecondsPerSweep, ...
    allRiseTimeInMilliSecondsPerSweep];


%% CELL ANALYSIS - summed PSC ==============================================

% sum all the ROI averages
summedCurrent = sum(meanDataInROI,2);

% adjusting pulseOnset value if this is a problem file
% remember that xmin = lightOnsetTime - baselineDurationInSeconds;
if problemFile == 1
    pulseOnset = baselineDurationInSeconds * samplingFrequency;
    afterLightDataPoint = pulseOnset + lightPulseAnalysisWindowInDataPts;
end    

% get the kinetics of the summed PSC - to compare with kinetics of PSC
% evoked by whole field illumination
summedCurrentOnsetLatency = min(find(conv2(sign(diff(summedCurrent(pulseOnset:afterLightDataPoint))), zeros(thresholdInDataPts,1)+(inwardORoutward), 'valid')==thresholdInDataPts));
% if onset was detected, get amplitude and location (index) of summedCurrent
% also calculate rise time
if ~isempty(summedCurrentOnsetLatency)      
    % look for peak current AFTER onset to avoid misleading riseTimes
    onsetDataPoint = pulseOnset + summedCurrentOnsetLatency;
    afterOnsetDataPoint = onsetDataPoint + lightPulseAnalysisWindowInDataPts;                
    if inwardORoutward == 1 
        [summedCurrentAmp, summedCurrentLoc] = max(summedCurrent(onsetDataPoint:afterOnsetDataPoint));
        summedCurrentTimeTo10percentOfPeakInDataPoints = find(summedCurrent(onsetDataPoint:afterOnsetDataPoint) >= 0.1 * summedCurrentAmp, 1);
        summedCurrentTimeTo90percentOfPeakInDataPoints = find(summedCurrent(onsetDataPoint:afterOnsetDataPoint) >= 0.9 * summedCurrentAmp, 1);
    else
        [summedCurrentAmp, summedCurrentLoc] = min(summedCurrent(onsetDataPoint:afterOnsetDataPoint));
        summedCurrentTimeTo10percentOfPeakInDataPoints = find(summedCurrent(onsetDataPoint:afterOnsetDataPoint) <= 0.1 * summedCurrentAmp, 1);
        summedCurrentTimeTo90percentOfPeakInDataPoints = find(summedCurrent(onsetDataPoint:afterOnsetDataPoint) <= 0.9 * summedCurrentAmp, 1);
    end         
% if onset was not detected, NaN it all! 
else
    summedCurrentOnsetLatency = NaN;
    summedCurrentAmp = NaN;
    summedCurrentLoc = NaN;
    summedCurrentTimeTo10percentOfPeakInDataPoints = NaN;
    summedCurrentTimeTo90percentOfPeakInDataPoints = NaN;
end

% Convert data points to milliseconds (1000 multiplication is conversion from seconds to milliseconds)
summedCurrentOnsetLatencyInMilliSeconds = 1000*summedCurrentOnsetLatency/samplingFrequency;
summedCurrentPeakLatencyInMilliSeconds = 1000*summedCurrentLoc/samplingFrequency;
summedCurrentTimeTo10percentOfPeakInMilliSeconds = 1000*summedCurrentTimeTo10percentOfPeakInDataPoints/samplingFrequency;
summedCurrentTimeTo90percentOfPeakInMilliSeconds = 1000*summedCurrentTimeTo90percentOfPeakInDataPoints/samplingFrequency;
summedCurrentRiseTimeInMilliSeconds = summedCurrentTimeTo90percentOfPeakInMilliSeconds - summedCurrentTimeTo10percentOfPeakInMilliSeconds; 


%% CELL ANALYSIS - charge per ROI
% average area under the curve according to Gordon Sheppard's suggestion

% meanDataInROI should have ROI columns.

size(meanDataInROI);

% gather current in response to light from light onset to light onset +
% user defined interval (15 ms if lightPulseAnalysisWindowInSeconds = 0.015)
dataSubsetForCharge = meanDataInROI(pulseOnset:afterLightDataPoint,:);
size(dataSubsetForCharge);

% find area under the curve using trapz
% use dim 2 to integrate over each column
avgChargeInROI = trapz(dataSubsetForCharge,1);

% adjust charge value according to sampling rate
avgChargeInROI = avgChargeInROI * 1/samplingFrequency;


%% CELL ANALYSIS - exclude "PerROI" data in ROIs with 100% failures ==============================================

% get success rate
successRatePerROI = (100 - failureRatePerROI)/100;

% find ROIs in which success rate is 0%
% aka according to our criteria, there are not light-evoked events in this ROI
failedROIs = find(successRatePerROI==0);

% remove data from other "PerROI" variables 
% this will yield cleaner heatmaps
meanPeakPerROI_noFailures = meanPeakPerROI;
meanPeakPerROI_noFailures(failedROIs) = NaN;

meanOnsetLatencyPerROI_noFailures = meanOnsetLatencyPerROI;
meanOnsetLatencyPerROI_noFailures(failedROIs) = NaN;

meanRiseTimePerROI_noFailures = meanRiseTimePerROI;
meanRiseTimePerROI_noFailures(failedROIs) = NaN;

stdOnsetLatencyPerROI_noFailures = stdOnsetLatencyPerROI;
stdOnsetLatencyPerROI_noFailures(failedROIs) = NaN;


%% CELL ANALYSIS - exclude "PerROI" data in ROIs with >50% failures ==============================================
% this is akin to one of the criteria in Ambrosi et al for classifying a
% cells as "shows an oIPSC" - here's the quote from the methods: "In rare sweeps,
% mIPSCs were mislabeled as oIPSCs. Thus, a cell was labeled as ‘‘shows an
% oIPSC’’ only if oIPSCs were detected in more than 50% of the recorded
% sweeps. Cells that did not fit this criteria were labeled as ‘‘no
% oIPSC’."

if discardROIsWithLowFreq == 1
    % find ROIs in which success rate is <50%
    % aka according to this criteria, there are not light-evoked events in this ROI
    % ALERT: might want to soft code the failure threshold! Instead of
    % discarding ROIs with > 50% failures, I might want to change that to
    % >75% failures or something, as I increase the number of sweeps per
    % ROI!
    failedROIs = find(successRatePerROI<0.5);

    % remove data from "PerROI" variables 
    meanPeakPerROI_noFailures(failedROIs) = NaN;
    meanOnsetLatencyPerROI_noFailures(failedROIs) = NaN;
    meanRiseTimePerROI_noFailures(failedROIs) = NaN;
    stdOnsetLatencyPerROI_noFailures(failedROIs) = NaN;    
end


%% CELL ANALYSIS - data storage =================================================================

% store subset of analyzed sweeps
allAnalyzedSweeps = setdiff(allSweeps,discardedSweeps);

% Store cell-specific data
dataCell = [mouseNumber, ...
    experimentDate, ...
    firstSweepNumber, ...
    lastSweepNumber, ...    
    gridColumns, ...
    gridRows, ...
    orderedGrid, ...
    plannedSweepsPerROI, ...   
    length(discardedSweeps), ...
    length(allAnalyzedSweeps), ...
    lightChannel, ...
    ledPowerChannel, ...
    singleLightPulse, ...        
    inwardORoutward, ...
    baselineDurationInSeconds, ...
    lightPulseAnalysisWindowInSeconds, ...
    thresholdInDataPts, ...
    amplitudeThreshold, ...   
    smoothSpan, ...    
    discardROIsWithLowFreq, ...
    problemFile, ...
    rsTestPulseOnsetTime, ...
    autoRsOnsetTime, ...
    voltageCmdChannel, ...
    stimDur, ... 
    stimFreq, ...
    lightDur, ...
    ledPower, ...
    mean(allRs, 'omitnan'), ...
    min(allRs), ...
    max(allRs), ...
    mean(baselineCurrentAll, 'omitnan'), ...
    min(baselineCurrentAll), ...
    max(baselineCurrentAll), ...
    summedCurrentAmp, ...
    summedCurrentOnsetLatencyInMilliSeconds, ...
    summedCurrentPeakLatencyInMilliSeconds, ...
    summedCurrentRiseTimeInMilliSeconds];

% store ROI-by-ROI data in multiple rows
dataCellMultipleRows = [repmat(dataCell,totalROIs,1), ...
    [1:totalROIs]', ...
    sweepsPerROI, ...
    successRatePerROI', ...
    meanPeakPerROI_noFailures', ...
    meanOnsetLatencyPerROI_noFailures', ...
    meanRiseTimePerROI_noFailures', ...
    stdOnsetLatencyPerROI_noFailures', ...
    avgChargeInROI'];

% note: I can't export sweepsInROI easily because ROIs might have different
% total number of sweeps

% % store ROI-by-ROI data in a single row
% dataCellSingleRow = dataCell;


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


%% PLOT 3 - Rs ===============================================================================================

% figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_batch - Rs')); % naming figure file
figure('name', strcat(fileName, 'rs')); % naming figure file
plot(allAnalyzedSweeps, allRs,'-o');
% plot lines marking 30% increase and 30% decrese in Rs compared to first test pulse
line([allAnalyzedSweeps(1) allAnalyzedSweeps(end)],[allRs(1)*0.7, allRs(1)*0.7],'Color','black','LineStyle','--')
line([allAnalyzedSweeps(1) allAnalyzedSweeps(end)],[allRs(1)*1.3, allRs(1)*1.3],'Color','black','LineStyle','--')
axis([allAnalyzedSweeps(1) inf 0 60])
ylabel('Rs (M\Omega)');
xlabel('Sweeps');
title([obj.file ' rs'],'Interpreter','none');
movegui('southeast');


%% PLOT 4 - Tiled oIPSC amplitude =========================================================================

% create figure & name it
figure('name', strcat(fileName, 'niceplot'));
t = tiledlayout(gridRows, gridColumns);

% plotting niceplots 
for ROI = 1:totalROIs
    nexttile
    hold on;
    
    % plot individual sweeps
    for sweep = cell2mat(relativeSweepsInROI(ROI))
        % if this is a problem file (aka the light onset is not consistent
        % accross sweeps), plot the data subset around the light pulse
        if problemFile ==1
            plot(xAroundLightPulse, yAroundLightPulseAll(:, sweep),'Color',[0, 0, 0, 0.25]);
        else
            plot(x, yBaselineSubAll(:, sweep),'Color',[0, 0, 0, 0.25]);
        end
    end
    
    % plot average acrross sweeps
%     plot(x, yBaselineSubAllMeanBySquare(:, square),'Color','black','LineWidth',0.7); 
%     line([xmin xmax],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--');
    axis([xmin xmax ymin ymax]);    
    
    % commented out on 2024-10-29 to prep paper figs
    % % adding light stim - individual pulses   
    % % ALERT: note that this code will use light stim parameters from the last sweep!
    % % if light stim is not the same accross all sweeps, this will be
    % % misleading!
    % for nStim=1:length(lightPulseStart)
    %     line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-inwardORoutward*(ymax),-inwardORoutward*(ymax)],'Color',[0 0.4470 0.7410],'LineWidth',5)
    % end
            
    % remove x and y labels from all ROIs
    xticklabels([]);
    yticklabels([]);
    
    % add scale bar to last plot
    if ROI == gridRows * gridColumns
        xmaxScale = xmax;
        xminScale = xmin;
        line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
        line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
%         text(xmaxScale-(xmaxScale-xminScale)/9,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
%         text(xmaxScale-(xmaxScale-xminScale)/7,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/7)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
        text(xmaxScale-(xmaxScale-xminScale),ymin+((ymax-ymin)/10),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
        text(xmaxScale-(xmaxScale-xminScale),ymin+((ymax-ymin)/3),strcat(num2str((ymax-ymin)/7)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))

    end
    
    hold off;    
    set(gca,'Visible','off');

end

t.TileSpacing = 'compact';
t.Padding = 'compact';
% xlabel(t,'Time (s)')
% ylabel(t,'oIPSC (pA)');
% ylabel(t, strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));

% set figure size to the same as the cropped cell image
% FYI this only adjusts the outer figure size, not the inner figure size...
% set(gcf,'OuterPosition',[1 1 outerWidth outerHeight]);
% set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);


%% PLOT 4.1 - Tiled oIPSC COLORS =========================================================================

% color blind diverging color scheme
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

% create figure & name it
figure('name', strcat(fileName, 'niceplot order'));
t = tiledlayout(gridRows, gridColumns);

% plotting niceplots 
for ROI = 1:totalROIs
    nexttile
    hold on;

    % had to add this to make the color plot work
    sweepNumFrom1to10 = 1;

    % plot individual sweeps
    for sweep = cell2mat(relativeSweepsInROI(ROI))
        % if this is a problem file (aka the light onset is not consistent
        % accross sweeps), plot the data subset around the light pulse
        if problemFile ==1
            plot(xAroundLightPulse, yAroundLightPulseAll(:, sweep),'Color',brownToPurpleRgb(sweepNumFrom1to10,:));
        else
            plot(x, yBaselineSubAll(:, sweep),'Color',brownToPurpleRgb(sweepNumFrom1to10,:));
        end
        
        % had to add this to make the color plot work
        sweepNumFrom1to10 = sweepNumFrom1to10 + 1;
    end
    
    % plot average acrross sweeps
%     plot(x, yBaselineSubAllMeanBySquare(:, square),'Color','black','LineWidth',0.7); 
%     line([xmin xmax],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--');
    axis([xmin xmax ymin ymax]);    
    
    % commented out on 2024-10-29 to prep paper figs
    % % adding light stim - individual pulses   
    % % ALERT: note that this code will use light stim parameters from the last sweep!
    % % if light stim is not the same accross all sweeps, this will be
    % % misleading!
    % for nStim=1:length(lightPulseStart)
    %     line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-inwardORoutward*(ymax),-inwardORoutward*(ymax)],'Color',[0 0.4470 0.7410],'LineWidth',5)
    % end
    
    % remove x and y labels from all ROIs
    xticklabels([]);
    yticklabels([]);
    
    % add scale bar to last plot
    if ROI == gridRows * gridColumns
        xmaxScale = xmax;
        xminScale = xmin;
        line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
        line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
%         text(xmaxScale-(xmaxScale-xminScale)/9,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
%         text(xmaxScale-(xmaxScale-xminScale)/7,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/7)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
        text(xmaxScale-(xmaxScale-xminScale),ymin+((ymax-ymin)/10),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
        text(xmaxScale-(xmaxScale-xminScale),ymin+((ymax-ymin)/3),strcat(num2str((ymax-ymin)/7)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))

        % add color legend?
        colormap(brownToPurpleRgb(1:max(sweepsPerROI),:))
        c = colorbar('Ticks',[1,max(sweepsPerROI)]);
        caxis([1 max(sweepsPerROI)]);
        c.Label.String = 'Sweep #';
        set(c, 'YDir', 'reverse');

    end
    
    hold off;    
    set(gca,'Visible','off');

end

t.TileSpacing = 'compact';
t.Padding = 'compact';
% xlabel(t,'Time (s)')
% ylabel(t,'oIPSC (pA)');
% ylabel(t, strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));

% set figure size to the same as the cropped cell image
% FYI this only adjusts the outer figure size, not the inner figure size...
% set(gcf,'OuterPosition',[1 1 outerWidth outerHeight]);
% set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);


%% PLOT 4.2 - tiled averages

% create figure & name it
figure('name', strcat(fileName, 'niceplot avgs'));
t = tiledlayout(gridRows, gridColumns);

% plotting niceplots 
for ROI = 1:totalROIs
    nexttile
    hold on;
    
    % plot avg
    if problemFile ==1
        plot(xAroundLightPulse, meanDataInROI(:, ROI),'Color',[0, 0, 0, 1]);
    else
        plot(x, meanDataInROI(:,ROI),'Color',[0, 0, 0, 1]);
    end
    
    axis([xmin xmax ymin ymax]);    
    
    % commented out on 2024-10-29 to prep paper figs
    % % adding light stim - individual pulses   
    % % ALERT: note that this code will use light stim parameters from the last sweep!
    % % if light stim is not the same accross all sweeps, this will be
    % % misleading!
    % for nStim=1:length(lightPulseStart)
    %     line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-inwardORoutward*(ymax),-inwardORoutward*(ymax)],'Color',[0 0.4470 0.7410],'LineWidth',5)
    % end
    
    % remove x and y labels from all ROIs
    xticklabels([]);
    yticklabels([]);
    
    % add scale bar to last plot
    if ROI == gridRows * gridColumns
        xmaxScale = xmax;
        xminScale = xmin;
        line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
        line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
        text(xmaxScale-(xmaxScale-xminScale),ymin+((ymax-ymin)/10),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
        text(xmaxScale-(xmaxScale-xminScale),ymin+((ymax-ymin)/3),strcat(num2str((ymax-ymin)/7)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))

    end
    
    hold off;    
    set(gca,'Visible','off');

end

t.TileSpacing = 'compact';
t.Padding = 'compact';

% set figure size to the same as the cropped cell image
% FYI this only adjusts the outer figure size, not the inner figure size...
% set(gcf,'OuterPosition',[1 1 outerWidth outerHeight]);
% set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);



%% PLOT 5 - subtracted baseline current =====================================================

figure('name', strcat(fileName, 'baseline current')); % naming figure file
plot(allAnalyzedSweeps, baselineCurrentAll,'-o');
axis([allAnalyzedSweeps(1) inf -500 500])
ylabel('Subtracted Baseline Current (pA)');
xlabel('Sweeps');
title([obj.file ' baseline current'],'Interpreter','none');
movegui('southwest');


%% PLOT 6 - summed PSCs

% match the xmin and xmax from psc_vs_light_single
xminHere = lightOnsetTime-0.02;
xmaxHere = lightOnsetTime+0.2;

figure('name', strcat(fileName, 'summed PSCs')); % naming figure file
hold on;
if problemFile ==1
    plot(xAroundLightPulse, meanDataInROI,'Color',[0, 0, 0, 0.25]);
    plot(xAroundLightPulse, summedCurrent,'Color','black','LineWidth',0.7); 
else
    plot(x, meanDataInROI,'Color',[0, 0, 0, 0.25]);
    plot(x, summedCurrent,'Color','black','LineWidth',0.7); 
end
line([xminHere xmaxHere],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--');
axis([xminHere xmaxHere ymin ymax]);

% adding light stim - individual pulses   
% note that this code will use light stim parameters from the last sweep!
% if light stim is not the same accross all sweeps, this will be
% misleading!
for nStim=1:length(lightPulseStart)
    line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-inwardORoutward*(ymax),-inwardORoutward*(ymax)],'Color',[0 0.4470 0.7410],'LineWidth',10)
end

% add scale bar
xmaxScale = xmaxHere;
xminScale = xminHere;
line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
text(xmaxScale-(xmaxScale-xminScale)/9,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
text(xmaxScale-(xmaxScale-xminScale)/7,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/7)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))

hold off;    
set(gca,'Visible','off');
set(gcf,'Position',[1400 550 500 400]);
movegui('south');


%% PLOT 7 - heatmap of success rate (0 to 100%)

% organize data for heatmap
dataForHeatmap = reshape(successRatePerROI,gridColumns,[]).';

% resize heatmap
resizedHeatmap = imresize(dataForHeatmap, [cropHeight cropWidth], 'nearest');

% make heatmap without the heatmap function
fig1 = figure('name', strcat(fileName, 'probability')); % naming figure file
ax1 = axes('Parent',fig1);
% imshow(resizedHeatmap,'Colormap',flipud(parula),'DisplayRange', [0,100], 'Border', 'tight');
imshow(resizedHeatmap,'Colormap',customColorMapPink,'DisplayRange', [0,1], 'Border', 'tight', 'Parent', ax1);
title(ax1, 'P(oIPSC)');
c1 = colorbar(ax1);
c1.Label.String = 'P(oIPSC)';


%% PLOT 9.2 - heatmap of average charge (AUC) normalized to largest charge

% set hetmap edges
heatmapMin = 0;
heatmapMax = 1;

% normalize meanPeakPerROI by largest peak
if inwardORoutward == 1
    maxCharge = max(avgChargeInROI);
else
    maxCharge = min(avgChargeInROI);
end
normalizedToMaxCharge = avgChargeInROI/maxCharge;

% organize data for heatmap
dataForHeatmap = reshape(normalizedToMaxCharge,gridColumns,[]).';

% resize heatmap
resizedHeatmap = imresize(dataForHeatmap, [cropHeight cropWidth], 'nearest');

% make heatmap without the heatmap function
fig2 = figure('name', strcat(fileName, 'charge')); % naming figure file
ax2 = axes('Parent',fig2);
imshow(resizedHeatmap,'Colormap',customColorMapVermillion,'DisplayRange', [heatmapMin,heatmapMax], 'Border', 'tight','Parent',ax2);
title(ax2, strcat('Normalized avg charge - max charge = ', num2str(round(maxCharge))));
c2 = colorbar(ax2);
c2.Label.String = 'Normalized AVG charge';


%% PLOT 10 - heatmap of onset latencies 

% organize data for heatmap
dataForHeatmap = reshape(meanOnsetLatencyPerROI_noFailures,gridColumns,[]).';

% resize heatmap
resizedHeatmap = imresize(dataForHeatmap, [cropHeight cropWidth], 'nearest');

% make heatmap without the heatmap function
fig3 = figure('name', strcat(fileName, 'latency')); % naming figure file
ax3 = axes('Parent',fig3);
% fig = imshow(resizedHeatmap,'Colormap',flipud(parula),'DisplayRange', [0,10], 'Border', 'tight');
fig3_1 = imshow(resizedHeatmap,'Colormap',flipud(customColorMapBlue),'DisplayRange', [0,10], 'Border', 'tight','Parent',ax3);
title(ax3,'oIPSC onset latency (ms)')
% set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
set(fig3_1, 'AlphaData', ~isnan(resizedHeatmap));
c3 = colorbar(ax3);
c3.Label.String = 'Onset Latency (ms)';


%% PLOT 11 - heatmap of onset latency jitter 

% organize data for heatmap
dataForHeatmap = reshape(stdOnsetLatencyPerROI_noFailures,gridColumns,[]).';

% resize heatmap
resizedHeatmap = imresize(dataForHeatmap, [cropHeight cropWidth], 'nearest');

% make heatmap without the heatmap function
fig4=figure('name', strcat(fileName, 'jitter')); % naming figure file
ax4 = axes('Parent',fig4);
% fig = imshow(resizedHeatmap,'Colormap',flipud(parula),'DisplayRange', [0,10], 'Border', 'tight');
fig4_1 = imshow(resizedHeatmap,'Colormap',flipud(customColorMapSky),'DisplayRange', [0,5], 'Border', 'tight','Parent',ax4);
title(ax4,'oIPSC onset latency jitter (ms)')
% set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
set(fig4_1, 'AlphaData', ~isnan(resizedHeatmap));
c4 = colorbar(ax4);
c4.Label.String = 'Onset Latency Jitter (ms)';


%% PLOT 12 - heatmap of rise time 

% organize data for heatmap
dataForHeatmap = reshape(meanRiseTimePerROI_noFailures,gridColumns,[]).';

% resize heatmap
resizedHeatmap = imresize(dataForHeatmap, [cropHeight cropWidth], 'nearest');

% make heatmap without the heatmap function
fig5=figure('name', strcat(fileName, 'rise')); % naming figure file
ax5=axes('Parent',fig5);
% fig = imshow(resizedHeatmap,'Colormap',flipud(parula),'DisplayRange', [0,10], 'Border', 'tight');
fig5_1 = imshow(resizedHeatmap,'Colormap',flipud(customColorMapSky),'DisplayRange', [0,5], 'Border', 'tight','Parent',ax5);
title(ax5,'oIPSC rise time (ms)')
% set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
set(fig5_1, 'AlphaData', ~isnan(resizedHeatmap));
c5 = colorbar(ax5);
c5.Label.String = 'Rise time (ms)';


%% EXPORTING XLS files - sweep-by-sweep ==========================================

% create cell array with strings for naming the amplitude and latency of the oPSC for each light pulse
oPSCvariableNamesCharge = cell(1,length(allChargesPerLightPulse));
oPSCvariableNamesFailure = cell(1,length(allFailureAssessments));
oPSCvariableNamesAmplitude = cell(1,length(allLightEvokedCurrentsAmp));
oPSCvariableNamesOnsetLatency = cell(1,length(allLightEvokedResponseOnsetLatencyInMilliSeconds));
oPSCvariableNamesPeakLatency = cell(1,length(allLightEvokedResponsePeakLatencyInMilliSeconds));
oPSCvariableNames10percentLatency = cell(1,length(allTimeTo10percentOfPeakInMilliSeconds));
oPSCvariableNames90percentLatency = cell(1,length(allTimeTo90percentOfPeakInMilliSeconds));
oPSCvariableNamesRiseTime = cell(1,length(allRiseTimeInMilliSeconds));
for lightPulse = 1:length(allLightEvokedCurrentsAmp)
    oPSCvariableNamesCharge(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(pA*s)')};
    oPSCvariableNamesFailure(lightPulse) = {strcat(num2str(lightPulse), 'isFail')};
    oPSCvariableNamesAmplitude(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(pA)')};
    oPSCvariableNamesOnsetLatency(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(onsetLat_ms)')};
    oPSCvariableNamesPeakLatency(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(peakLat_ms)')};
    oPSCvariableNames10percentLatency(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(10peakLat_ms)')};
    oPSCvariableNames90percentLatency(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(90peakLat_ms)')};
    oPSCvariableNamesRiseTime(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(riseTime_ms)')};
end

% store key sweep-by-sweep data
filename = strcat(fileName, 'sweep_by_sweep'); 
fulldirectory = strcat(filename,'.xls'); 
fulldirectory = fullfile(savefileto,fulldirectory); % fullfile adds OS-appropriate slash
dataInCellFormat = {};
dataInCellFormat = num2cell(data);
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
    {'mouse', ...
    'date', ...
    'sweep', ...
    'ROI', ...
    'gridColumns', ...
    'gridRows', ...
    'isOrderedGrid(1)orNot(0)', ...
    'plannedSweepsPerROI', ...
    'lightChannel', ...
    'ledPowerChannel', ...
    'singleLightPulse(1)orTrain(0)', ...
    'inward(-1)ORoutward(1)', ...
    'baselineDurationInSeconds', ...
    'lightPulseAnalysisWindowInSeconds', ...
    'thresholdInDataPts', ...
    'amplitudeThreshold(pA)', ...  
    'smoothSpan(dataPts)', ...
    'isProblemFile(1)orNot(0)', ...
    'rsTestPulseOnsetTime', ...
    'autoRsOnsetTime(1)orManual(0)', ...
    'voltageCmdChannel', ...
    'lightPulseDur(s)', ... 
    'lightStimFreq(Hz)', ...
    'lightDur(s)', ...   
    'ledPower(V)', ...   
    'seriesResistance(Mohm)', ...
    'baselineCurrent(pA)', ...
    'heavyEditing(1)orNot(0)', ...
    oPSCvariableNamesCharge{:}, ...
    oPSCvariableNamesFailure{:}, ...
    oPSCvariableNamesAmplitude{:}, ...
    oPSCvariableNamesOnsetLatency{:}, ...
    oPSCvariableNamesPeakLatency{:}, ...
    oPSCvariableNames10percentLatency{:}, ... 
    oPSCvariableNames90percentLatency{:}, ...
    oPSCvariableNamesRiseTime{:}});

writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the sweep_by_sweep xls file')


%% EXPORTING XLS files - ROI-by-ROI ==========================================

% ALERT: warn user about current code limitations
% cuz I don't want to deal with a gigantic csv containing multiple light
% pulses for multiple ROIs!
if singleLightPulse == 0
    disp('only analysing the first oIPSC in each sweep bro, look out if you got more')
end    

% % % create cell array for variable naming (not currently used because I can't
% % % easily store sweepsInROI
% % headings = cell(1, plannedSweepsPerROI);
% % for sweep = [1:plannedSweepsPerROI]
% %     headings(sweep) = {strcat('sweep', num2str(sweep))};
% % end

% store key ROI-by-ROI data - each row is a square
filename = strcat(fileName, 'ROI_by_ROI'); 
fulldirectory = strcat(filename,'.xls'); 
fulldirectory = fullfile(savefileto,fulldirectory); % fullfile adds OS-appropriate slash
dataInCellFormat = {};
dataInCellFormat = num2cell(dataCellMultipleRows);
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
    {'mouse', ...
    'date', ...
    'firstSweep', ...
    'lastSweep', ...
    'gridColumns', ...
    'gridRows', ...
    'isOrderedGrid(1)orNot(0)', ...
    'plannedSweepsPerROI', ...   
    'nDiscardedSweeps', ...
    'nAnalyzedSweeps', ...
    'lightChannel', ...
    'ledChannel', ...
    'singleLightPulse(1)orTrain(0)', ...
    'inward(-1)ORoutward(1)', ...
    'baselineDurationInSeconds', ...
    'lightPulseAnalysisWindowInSeconds', ...
    'thresholdInDataPts', ...
    'amplitudeThreshold(pA)', ...
    'smoothSpan(dataPts)', ...
    'discardROIsWithLowFreq', ...
    'isProblemFile(1)orNot(0)', ...
    'rsTestPulseOnsetTime', ...
    'autoRsOnsetTime(1)orManual(0)', ...
    'voltageCmdChannel', ...
    'lightPulseDur(s)', ... 
    'lightStimFreq(Hz)', ...
    'lightDur(s)', ...   
    'ledPower(V)', ... 
    'AVGseriesResistance(Mohm)', ...
    'MINseriesResistance(Mohm)', ...
    'MAXseriesResistance(Mohm)', ...
    'AVGbaselineCurrent(pA)', ...
    'MINbaselineCurrent(pA)', ...
    'MAXbaselineCurrent(pA)', ...
    'summedCurrentAmp(pA)', ...
    'summedCurrentOnsetLatency(ms)', ...
    'summedCurrentPeakLatency(ms)', ...
    'summedCurrentRiseTime(ms)', ...
    'ROI', ...
    'sweepsPerROI', ...
    'successRate', ...
    'AVGpeakAmp(pA)', ...
    'AVGonsetLatency(ms)', ...
    'AVGriseTime(ms)', ...
    'SDonsetLatency(ms)', ...
    'avgChargeInROI(pA*ms)'});
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the ROI_by_ROI xls file')

end
