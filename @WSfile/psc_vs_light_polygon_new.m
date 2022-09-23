%{ 
DOCUMENTATION
Created: 2022 09 21
Works? Yes
Author: PA

This function is used to analyze light-evoked PSCs in sCRACM experiments using the polygon.
This function was derived from psc_vs_light_polygon, but went through A LOT of rewritting.
The data organization is very different, to better handle discarded sweeps
without losing track of which sweeps belong to which ROIs (Regions Of
Interest - subarea of the field that was illuminated).

INPUTS explained:
    - gridColumns: number of columns in the polygon ROI grid

    - gridRows: number of rows in the polygon ROI grid

    - discardedSweeps: specific sweeps that should NOT be analyzed due to
    artifacts. If I want to discard sweep 0024, just write 24.

    - lightChannel: channel where the info about the light stim is stored.
    Usually 2 (LED 1), 3 (LED 2), or 4 (polygon). 

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
    about sign (+ or -). The default 25 chosen based on the noise level of
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

    - savefileto: save CSV files to this folder.

INPUTS defaults:
    % Affects data analysis - Organizing data by o-stim grid
    gridColumns = 5;
    gridRows = 5;

    % Affects data analysis - Finding/quantifyting oIPSCs
    discardedSweeps = [];
    discardedSweepsFromEnd = 0;
    lightChannel = 4;
    singleLightPulse = 1; 
    inwardORoutward = -1;               
    baselineDurationInSeconds = 0.01;
    lightPulseAnalysisWindowInSeconds = 0.02;
    thresholdInDataPts = 5;     
    amplitudeThreshold = 25;
    smoothSpan = 3; 

    % Affects data analysis - Calculating Rs
    rsTestPulseOnsetTime = 1;
    autoRsOnsetTime = 0;
    voltageCmdChannel = 2;

    % Affects data display: 
    ymin = -1500;           %-2050
    ymax = 600;             %50
    cellImageFileNameDIC = 's1c1_z1_moved1_dic.tif';
    cellImageFileNameAlexa = 's1c1_z1_moved1_647.tif';
    cellImageDir = 'D:\NU server\Priscilla - BACKUP 20200319\Ephys\2022\20220728 m076 dms loop';

    % Affects data saving:
    savefileto = 'D:\Temp\From MATLAB 2022 08 31 psc';    
    
OUTPUTS:
    Fig
    XLS

ASSUMPTIONS: 
    - An ordered polygon grid design was used - aka not random
    - Width and frequency of polygon stim is the the same as LED o-stim
    (only relevant if you use the polygon channel as the light stim
    channel)

BEWARE:
    - ctrl+f for ALERT and ASSUMPTION

TO DO:
    - noticed weird bug on 8/31/22: if I include the last sweep on the
    discarded sweeps list, matlab does not plot any of the tiled niceplots!
%}

% obj = m729.s0246;
% obj = m729.s0111;
% obj = m729.s0371;
% obj = m731.s0007

function psc_vs_light_polygon_new(obj)
%%  USER INPUT ============================================================

% Affects data analysis - Organizing data by o-stim grid
gridColumns = 7;
gridRows = 7;

% Affects data analysis - Finding/quantifyting oIPSCs
discardedSweeps = [];
lightChannel = 4;
ledPowerChannel = 3;
singleLightPulse = 1; 
inwardORoutward = -1;                       % 1 (positive) is outward; -1 (negative) in inward
baselineDurationInSeconds = 0.01;
lightPulseAnalysisWindowInSeconds = 0.015;  % ALERT: changed from 0.02 to 0.01 to 0.015 on 2022-09-21
thresholdInDataPts = 10;                    % ALERT! Changed from 10 to 5 to 10 (2022-09-23)
amplitudeThreshold = 25;                    % ALERT! this is new (2022-09-21)
smoothSpan = 3;                             % ALERT! this is new (2022-09-23)

% Affects data analysis - Calculating Rs
rsTestPulseOnsetTime = 1;
autoRsOnsetTime = 0;
voltageCmdChannel = 2;

% Affects data display: 
ymin = -3600;           %-2050      -3600
ymax = 600;             %50         600
cellImageFileNameDIC = 's2c4_z1_dic.tif';
cellImageFileNameAlexa = 's2c4_MAX_Stack Rendered Paths.tif';
cellImageDir = 'D:\NU server\Priscilla - BACKUP 20200319\Ephys\2022\20220914 m729 asc spiral';

% Affects data saving:
savefileto = 'D:\Temp\From MATLAB 2022 09 23 filtered';


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

% calculating variables based on user input - the Rs-related variables
% might change later in the code depending on the value in autoRsOnsetTime
baselineDurationInDataPts = baselineDurationInSeconds * samplingFrequency;
lightPulseAnalysisWindowInDataPts = lightPulseAnalysisWindowInSeconds * samplingFrequency;
rsBaselineDataPointInterval = ((rsTestPulseOnsetTime-0.05)*samplingFrequency):(rsTestPulseOnsetTime*samplingFrequency);
rsFirstTransientDataPointInterval = (rsTestPulseOnsetTime*samplingFrequency):(rsTestPulseOnsetTime+0.0025)*samplingFrequency;

% add 'subset' to fileName in case discardedSweeps is not empty
% to prevent overwritting of saved files with the whole dataset
if ~isempty(discardedSweeps)
    fileName = strcat(obj.file, '_subset');
end

% getting sweep numbers from file name
[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj); 

% ALERT: NEED TO TEST
% checking for incomplete sweeps and adding them to the list of discarded sweeps
% to avoid this error: "Index exceeds the number of array elements (0)". 
if numel(fieldnames(obj.sweeps)) < obj.header.NSweepsPerRun  
    lastCompleteSweep = firstSweepNumber + numel(fieldnames(obj.sweeps)) - 1;
    discardedSweeps = [discardedSweeps, lastCompleteSweep:lastSweepNumber];
end

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

% creating matrixes/arrays that will be filled later
allRs = [];
yBaselineSubAll = [];
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

% clearing cell arrays for safety
% (if you use this code as a script, data will not linger from one analysis
% to the next)
meanDataInROI = [];
peaksInROI = {};
onsetLatenciesInROI = {};
riseTimesInROI = {};
failuresInROI = {};
meanPeakPerROI = [];
meanOnsetLatencyPerROI = [];
meanRiseTimePerROI = [];
failureRatePerROI = [];
stdOnsetLatencyPerROI = [];


%% ROI BY ROI AND SWEEP BY SWEEP ANALYSIS ===================================================================================

% get data from sweeps in file (only the subset we will analyze)
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

        % this code that makes more sense
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
        
        % smooth data with moving average filter
        y = smooth(y,smoothSpan);

        % baseline subtraction
        baselineStart = lightPulseStart(1) - baselineDurationInDataPts;
        yBaselineSub = y-mean(y(baselineStart:lightPulseStart(1)));
        baselineCurrent = mean(y(baselineStart:lightPulseStart(1)));
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
            lightEvokedResponseOnsetLatencyInDataPoints = min(find(conv2(sign(diff(y(pulseOnset:afterLightDataPoint))), zeros(thresholdInDataPts,1)+(inwardORoutward), 'valid')==thresholdInDataPts));
            
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
            % light-evoked response in this sweep or if we whould consider
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
        dataInROI(ROI,column) = {yBaselineSub};
        
        % store key data for statistics (mean, SD, median) later
        peaksInROI(ROI,column) = {allLightEvokedCurrentsAmp};
        onsetLatenciesInROI(ROI,column) = {allLightEvokedResponseOnsetLatencyInMilliSeconds};
        riseTimesInROI(ROI,column) = {allRiseTimeInMilliSeconds};
        failuresInROI(ROI,column) = {isFailure};
        
        % store key data that will be exported (each row is a sweep)
        data = [data; ...
            mouseNumber, ...
            experimentDate, ...
            sweepNumber, ...         
            lightChannel, ...
            ledPowerChannel, ...
            singleLightPulse, ...        
            inwardORoutward, ...
            baselineDurationInSeconds, ...
            lightPulseAnalysisWindowInSeconds, ...
            thresholdInDataPts, ...
            amplitudeThreshold, ...
            rsTestPulseOnsetTime, ...
            autoRsOnsetTime, ...
            voltageCmdChannel, ...
            stimDur, ... 
            stimFreq, ...
            lightDur, ...
            ledPower, ...
            gridColumns, ...
            gridRows, ...
            plannedSweepsPerROI, ...
            seriesResistance, ...
            baselineCurrent]; 
        
    end 
    
    meanDataInROI = [meanDataInROI, mean(cell2mat(dataInROI(ROI,:)),2)];
    
    % ALERT: the code below might break down if you have multiple light
    % pulses, so I'm throwing a message out there
    if singleLightPulse ~= 1
        disp('bruh the "perROI" stuff was not made to handle multiple light pulses')
    end

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


%% CELL ANALYSIS - data storage =================================================================

% 
allAnalyzedSweeps = setdiff(allSweeps,discardedSweeps);

% Store cell-specific data
dataCell = [mouseNumber, ...
    experimentDate, ...
    firstSweepNumber, ...
    lastSweepNumber, ...
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
    rsTestPulseOnsetTime, ...
    autoRsOnsetTime, ...
    voltageCmdChannel, ...
    stimDur, ... 
    stimFreq, ...
    lightDur, ...
    ledPower, ...
    gridColumns, ...
    gridRows, ...
    plannedSweepsPerROI, ...
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
    stdOnsetLatencyPerROI_noFailures'];

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


%% PLOT 1 - cropped DIC cell image ===============================================================================
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


%% PLOT 3 - cropped Alexa cell image with grid ===============================================================================
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


%% PLOT 5 - Rs ===============================================================================================

figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_polygon - Rs all')); % naming figure file
plot(allAnalyzedSweeps, allRs,'-o');
% plot lines marking 30% increase and 30% decrese in Rs compared to first test pulse
line([allAnalyzedSweeps(1) allAnalyzedSweeps(end)],[allRs(1)*0.7, allRs(1)*0.7],'Color','black','LineStyle','--')
line([allAnalyzedSweeps(1) allAnalyzedSweeps(end)],[allRs(1)*1.3, allRs(1)*1.3],'Color','black','LineStyle','--')
axis([allAnalyzedSweeps(1) inf 0 60])
ylabel('Rs (M\Omega)');
xlabel('Sweeps');
title([obj.file ' rs'],'Interpreter','none');
movegui('southeast');


%% PLOT 6 - Tiled oIPSC amplitude =========================================================================

% create figure & name it
figure('name', strcat(fileName, '_', analysisDate, ' - psc_vs_light_polygon - tiled niceplots'));
t = tiledlayout(gridRows, gridColumns);

% plotting histogram plot 
for ROI = 1:totalROIs
    nexttile
    hold on;
    
    % plot individual sweeps
    for sweep = cell2mat(relativeSweepsInROI(ROI))
        plot(x, yBaselineSubAll(:, sweep),'Color',[0, 0, 0, 0.25]);
    end
    
    % plot average acrross sweeps
%     plot(x, yBaselineSubAllMeanBySquare(:, square),'Color','black','LineWidth',0.7); 
%     line([xmin xmax],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--');
    axis([xmin xmax ymin ymax]);    
    
    % adding light stim - individual pulses   
    % note that this code will use light stim parameters from the last sweep!
    % if light stim is not the same accross all sweeps, this will be
    % misleading!
    for nStim=1:length(lightPulseStart)
        line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-inwardORoutward*(ymax),-inwardORoutward*(ymax)],'Color',[0 0.4470 0.7410],'LineWidth',5)
    end
        
%     % remove y labels from all plots except the first
%     if square ~= 1 
%         yticklabels([]);
%     end
% 
%     % remove x labels from all plots except the last
%     if square ~= gridRows * gridColumns
%         xticklabels([]);
%     end
    
    % remove x and y labels from all squares
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
set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);


%% PLOT 7 - subtracted baseline current =====================================================

figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_polygon - baseline current')); % naming figure file
plot(allAnalyzedSweeps, baselineCurrentAll,'-o');
axis([allAnalyzedSweeps(1) inf -500 500])
ylabel('Subtracted Baseline Current (pA)');
xlabel('Sweeps');
title([obj.file ' baseline current'],'Interpreter','none');
movegui('southwest');


%% PLOT 8 - summed PSCs

% match the xmin and xmax from psc_vs_light_single
xminHere = lightOnsetTime-0.02;
xmaxHere = lightOnsetTime+0.2;

figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_polygon - summed PSCs')); % naming figure file
hold on;
plot(x, meanDataInROI,'Color',[0, 0, 0, 0.25]);
plot(x, summedCurrent,'Color','black','LineWidth',0.7); 
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


%% PLOT 12 - heatmap of success rate

% organize data for heatmap
dataForHeatmap = reshape(successRatePerROI,gridColumns,[]).';

% resize heatmap
resizedHeatmap = imresize(dataForHeatmap, [size(croppedImage,1) size(croppedImage,2)], 'nearest');

% make heatmap without the heatmap function
figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_polygon - success heatmap')); % naming figure file
% imshow(resizedHeatmap,'Colormap',flipud(parula),'DisplayRange', [0,100], 'Border', 'tight');
imshow(resizedHeatmap,'Colormap',customColorMapPink,'DisplayRange', [0,1], 'Border', 'tight');
title('p(oIPSC)');
set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
c = colorbar;
c.Label.String = 'p(oIPSC)';


%% PLOT 10 - heatmap of normalized PSCs (normalized to largest PSC)

% set hetmap edges
heatmapMin = 0;
heatmapMax = 1;

% normalize average peak-by-square by largest peak
% ALERT might want to change calculation of maxPSC to use meanPeakPerROI
% instead of allLightEvokedCurrentsPerSweep - cuz the max meanPeakPerROI will most
% likely be lower than the max allLightEvokedCurrentsPerSweep - DONE
if inwardORoutward == 1
    maxPSC = max(meanPeakPerROI);
else
    maxPSC = min(meanPeakPerROI);
end
normalizedToMaxPSC = meanPeakPerROI/maxPSC;

% exclude data from ROIs with 0% success rate
normalizedToMaxPSC_noFailures = normalizedToMaxPSC;
normalizedToMaxPSC_noFailures(failedROIs) = NaN;

% organize data for heatmap
dataForHeatmap = reshape(normalizedToMaxPSC_noFailures,gridColumns,[]).';

% resize heatmap
resizedHeatmap = imresize(dataForHeatmap, [size(croppedImage,1) size(croppedImage,2)], 'nearest');

% make heatmap without the heatmap function
figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_polygon - amplitude heatmap')); % naming figure file
% imshow(resizedHeatmap,'Colormap',parula,'DisplayRange', [heatmapMin,heatmapMax], 'Border', 'tight');
imshow(resizedHeatmap,'Colormap',customColorMapVermillion,'DisplayRange', [heatmapMin,heatmapMax], 'Border', 'tight');
title('Normalized oIPSC amplitude')
set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
c = colorbar;
c.Label.String = 'Normalized oIPSC amplitude';


%% PLOT 11.1 - heatmap of onset latencies 

% organize data for heatmap
dataForHeatmap = reshape(meanOnsetLatencyPerROI_noFailures,gridColumns,[]).';

% resize heatmap
resizedHeatmap = imresize(dataForHeatmap, [size(croppedImage,1) size(croppedImage,2)], 'nearest');

% make heatmap without the heatmap function
figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_polygon - latency heatmap')); % naming figure file
% fig = imshow(resizedHeatmap,'Colormap',flipud(parula),'DisplayRange', [0,10], 'Border', 'tight');
fig = imshow(resizedHeatmap,'Colormap',flipud(customColorMapBlue),'DisplayRange', [0,10], 'Border', 'tight');
title('oIPSC onset latency')
set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
set(fig, 'AlphaData', ~isnan(resizedHeatmap));
c = colorbar;
c.Label.String = 'Onset Latency';


%% PLOT 11.1 - heatmap of onset latency jitter 

% organize data for heatmap
dataForHeatmap = reshape(stdOnsetLatencyPerROI_noFailures,gridColumns,[]).';

% resize heatmap
resizedHeatmap = imresize(dataForHeatmap, [size(croppedImage,1) size(croppedImage,2)], 'nearest');

% make heatmap without the heatmap function
figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_polygon - jitter heatmap')); % naming figure file
% fig = imshow(resizedHeatmap,'Colormap',flipud(parula),'DisplayRange', [0,10], 'Border', 'tight');
fig = imshow(resizedHeatmap,'Colormap',flipud(customColorMapSky),'DisplayRange', [0,5], 'Border', 'tight');
title('oIPSC onset latency jitter')
set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
set(fig, 'AlphaData', ~isnan(resizedHeatmap));
c = colorbar;
c.Label.String = 'Onset Latency Jitter';


%% PLOT 11.1 - heatmap of rise time 

% organize data for heatmap
dataForHeatmap = reshape(meanRiseTimePerROI_noFailures,gridColumns,[]).';

% resize heatmap
resizedHeatmap = imresize(dataForHeatmap, [size(croppedImage,1) size(croppedImage,2)], 'nearest');

% make heatmap without the heatmap function
figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_polygon - rise heatmap')); % naming figure file
% fig = imshow(resizedHeatmap,'Colormap',flipud(parula),'DisplayRange', [0,10], 'Border', 'tight');
fig = imshow(resizedHeatmap,'Colormap',flipud(customColorMapSky),'DisplayRange', [0,5], 'Border', 'tight');
title('oIPSC rise time')
set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
set(fig, 'AlphaData', ~isnan(resizedHeatmap));
c = colorbar;
c.Label.String = 'Rise time (ms)';


%% EXPORTING XLS files - sweep-by-sweep ==========================================

% create cell array with strings for naming the amplitude and latency of the oPSC for each light pulse
oPSCvariableNamesFailure = cell(1,length(allFailureAssessments));
oPSCvariableNamesAmplitude = cell(1,length(allLightEvokedCurrentsAmp));
oPSCvariableNamesOnsetLatency = cell(1,length(allLightEvokedResponseOnsetLatencyInMilliSeconds));
oPSCvariableNamesPeakLatency = cell(1,length(allLightEvokedResponsePeakLatencyInMilliSeconds));
oPSCvariableNames10percentLatency = cell(1,length(allTimeTo10percentOfPeakInMilliSeconds));
oPSCvariableNames90percentLatency = cell(1,length(allTimeTo90percentOfPeakInMilliSeconds));
oPSCvariableNamesRiseTime = cell(1,length(allRiseTimeInMilliSeconds));
for lightPulse = 1:length(allLightEvokedCurrentsAmp)
    oPSCvariableNamesFailure(lightPulse) = {strcat(num2str(lightPulse), 'isFail')};
    oPSCvariableNamesAmplitude(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(pA)')};
    oPSCvariableNamesOnsetLatency(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(onsetLat_ms)')};
    oPSCvariableNamesPeakLatency(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(peakLat_ms)')};
    oPSCvariableNames10percentLatency(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(10peakLat_ms)')};
    oPSCvariableNames90percentLatency(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(90peakLat_ms)')};
    oPSCvariableNamesRiseTime(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(riseTime_ms)')};
end

% store key sweep-by-sweep data
filename = strcat(fileName, '_', analysisDate, " - psc_vs_light_polygon - sweep_by_sweep");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(data);
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
    {'mouse', ...
    'date', ...
    'sweep', ...
    'lightChannel', ...
    'ledPowerChannel', ...
    'singleLightPulse(1)orTrain(0)', ...
    'inward(-1)ORoutward(1)', ...
    'baselineDurationInSeconds', ...
    'lightPulseAnalysisWindowInSeconds', ...
    'thresholdInDataPts', ...
    'amplitudeThreshold(pA)', ...  
    'rsTestPulseOnsetTime', ...
    'autoRsOnsetTime(1)orManual(0)', ...
    'voltageCmdChannel', ...
    'lightPulseDur(s)', ... 
    'lightStimFreq(Hz)', ...
    'lightDur(s)', ...   
    'ledPower(V)', ...
    'gridColumns', ...
    'gridRows', ...
    'plannedSweepsPerROI', ...      
    'seriesResistance(Mohm)', ...
    'baselineCurrent(pA)', ...
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
filename = strcat(fileName, '_', analysisDate, " - psc_vs_light_polygon - ROI_by_ROI");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(dataCellMultipleRows);
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
    {'mouse', ...
    'date', ...
    'firstSweep', ...
    'lastSweep', ...
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
    'rsTestPulseOnsetTime', ...
    'autoRsOnsetTime(1)orManual(0)', ...
    'voltageCmdChannel', ...
    'lightPulseDur(s)', ... 
    'lightStimFreq(Hz)', ...
    'lightDur(s)', ...   
    'ledPower(V)', ...
    'gridColumns', ...
    'gridRows', ...
    'plannedSweepsPerROI', ...    
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
    'SDonsetLatency(ms)'});
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the ROI_by_ROI xls file')

end


%% Legacy code - not currently used

% % %% PLOT 9 - heatmap of normalized PSCs (normalized to summed PSC)
% % 
% % % set hetmap edges
% % heatmapMin = 0;
% % heatmapMax = 1;
% % 
% % % peak current of summed PSCs
% % summedPeak = min(summedPSC(pulseOnset:afterLightDataPoint))
% % 
% % % normalize average peak-by-square by summed peak
% % normalizedPSC = min(yMeanPerROI(pulseOnset:afterLightDataPoint,:))/min(summedPSC(pulseOnset:afterLightDataPoint));
% % 
% % % organize data for heatmap
% % dataForHeatmap = reshape(normalizedPSC,gridColumns,[]).';
% % 
% % % resize heatmap
% % resizedHeatmap = imresize(dataForHeatmap, [size(croppedImage,1) size(croppedImage,2)], 'nearest');
% % 
% % % make heatmap without the heatmap function
% % figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_polygon - heatmap 1')); % naming figure file
% % imshow(resizedHeatmap,'Colormap',parula,'DisplayRange', [heatmapMin,heatmapMax], 'Border', 'tight');
% % set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
% % c = colorbar;
% % c.Label.String = 'Normalized PSC';


% % %% PLOT 11.2 - heatmap of normalized onset latencies 
% % % this is a really confusing plot - do not recommend!
% % 
% % % normalize to longest latency
% % maxLatency = max(meanOnsetLatencyPerROI);
% % normalizedToMaxLatency = meanOnsetLatencyPerROI/maxLatency;
% % 
% % % organize data for heatmap
% % dataForHeatmap = reshape(normalizedToMaxLatency,gridColumns,[]).';
% % 
% % % resize heatmap
% % resizedHeatmap = imresize(dataForHeatmap, [size(croppedImage,1) size(croppedImage,2)], 'nearest');
% % 
% % % make heatmap without the heatmap function
% % figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_polygon - latency heatmap')); % naming figure file
% % % fig = imshow(resizedHeatmap,'Colormap',flipud(parula),'DisplayRange', [0,10], 'Border', 'tight');
% % fig = imshow(resizedHeatmap,'Colormap',flipud(customColorMap),'DisplayRange', [0,1], 'Border', 'tight');
% % title('oIPSC onset latency')
% % set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
% % set(fig, 'AlphaData', ~isnan(resizedHeatmap));
% % c = colorbar;
% % c.Label.String = 'Onset Latency';


% % %% PLOT 12 - heatmap of failure rate
% % % discontinued cuz success rate is more informative and fits my color
% % % coding better
% % 
% % % organize data for heatmap
% % dataForHeatmap = reshape(failureRatePerROI,gridColumns,[]).';
% % 
% % % resize heatmap
% % resizedHeatmap = imresize(dataForHeatmap, [size(croppedImage,1) size(croppedImage,2)], 'nearest');
% % 
% % % make heatmap without the heatmap function
% % figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_polygon - failure heatmap')); % naming figure file
% % imshow(resizedHeatmap,'Colormap',flipud(parula),'DisplayRange', [0,100], 'Border', 'tight');
% % set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
% % c = colorbar;
% % c.Label.String = 'Failure Rate (%)';


% % %% PLOT 1 - cropped DIC cell image ===============================================================================
% % % discontinued cuz I like the ones with the grid better
% % % make sure you're getting the image taken with zoom = 1
% % % concatenate strings from user input to get full path to figure file
% % cellImageFileDir = strcat(cellImageDir,'\',cellImageFileNameDIC);
% % 
% % % import image
% % cellImage = imread(cellImageFileDir);
% % 
% % % crop image according to the polygon's mirror calibration
% % % I verified this in PPT lol
% % % the original image is 1376 x 1024 pixels
% % % ASSUMPTION ALERT: the calibration of the polygon will not change over time
% % % I need to crop the top 100 pixels and the bottom 51 pixels
% % % imcrop determines the rectangle to keep in the following form: [XMIN YMIN WIDTH HEIGHT]
% % % note that y increases from top to bottom, so ymin should match my
% % % required "top crop".
% % % I do not need to crop along the x axis, so xmin = 1 and width = 1376
% % % the height is 1024 - 100 - 51 = 873
% % croppedImage = imcrop(cellImage, [1,100,1376,873]);
% % figure('name', strcat(fileName, '_', analysisDate, ' - psc_vs_light_polygon - DIC image'));
% % hold on;
% % imshow(croppedImage, 'Border', 'tight');
% % 
% % % calculate parameters for scale bar
% % % ASSUMPTION ALERT: pixelsPerMicron corresponds to my usual 40x objective at 1x zoom
% % xmaxImg = size(croppedImage,2);    % in pixels, should be 1376
% % ymaxImg = size(croppedImage,1);    % in pixels, should be 874
% % scaleDistanceFromEdge = 50;     % in pixels
% % scaleBarSize = 50;              % in um
% % pixelsPerMicron = 873 / 222.2;  % 222.2 um in 873 pixels
% % scaleBarSizeInPixels = scaleBarSize * pixelsPerMicron;
% % 
% % % add scale bar to figure
% % % line([x x], [y y])
% % line([scaleDistanceFromEdge scaleDistanceFromEdge+scaleBarSizeInPixels],...
% %     [ymaxImg-scaleDistanceFromEdge ymaxImg-scaleDistanceFromEdge],...
% %     'LineWidth', 2, 'Color', 'k');
% % % text(x, y, string)
% % text(scaleDistanceFromEdge,...
% %     ymaxImg-2*scaleDistanceFromEdge,...
% %     strcat(num2str(scaleBarSize), " μm"),...
% %     'FontSize', 10);
% % hold off;
% % 
% % % get inner figure size and store half of those values
% % pos = get(gcf, 'InnerPosition'); %// gives x left, y bottom, width, height
% % innerWidth = pos(3)/2;
% % innerHeight = pos(4)/2;
% % 
% % % get outer figure size and store half of those values
% % pos = get(gcf, 'OuterPosition'); %// gives x left, y bottom, width, height
% % outerWidth = pos(3)/2;
% % outerHeight = pos(4)/2;
% % 
% % % set figure size to the stored values
% % set(gcf,'InnerPosition',[0 maxHeight-innerHeight innerWidth innerHeight]);


% % %% EXPORTING XLS files - single row ==========================================
% % 
% % % store cell data with square-by-square info in a single row
% % filename = strcat(fileName, '_', analysisDate, " - psc_vs_light_polygon - square_by_square cell");
% % fulldirectory = strcat(savefileto,'\',filename,'.xls');        
% % dataInCellFormat = {};
% % dataInCellFormat = num2cell(dataCellSingleRow);
% % dataInCellFormat = [dataInCellFormat, cellImageFileNameDIC, cellImageFileNameAlexa];
% % 
% % % create cell array with strings for naming square-by-square data
% % % calculate the number of unique variables per square
% % numberOfVariablesPerSquare = size([sweepsPerSquarePerSquare, ...
% %                             mean(dataSubsetForMeanAndSD, 'omitnan'), ...
% %                             std(dataSubsetForMeanAndSD, 'omitnan')], 2);
% % % make a cell array with 1 row and as many columns as the number of unique
% % % variables we have (totalSquares*numberOfVariablesPerSquare)
% % squareVariables = cell(1,totalROIs*numberOfVariablesPerSquare);
% % % name each variable with square number and correct label
% % for ROI = [1:totalROIs]
% %     startingIndex = ROI*numberOfVariablesPerSquare+1 - numberOfVariablesPerSquare;
% %     squareVariables(startingIndex) = {strcat('sq', num2str(ROI), '-sweepsPerSquare')};
% %     squareVariables(startingIndex+1) = {strcat('sq', num2str(ROI), '-lightEvokedPeak(pA)AVG')};
% %     squareVariables(startingIndex+2) = {strcat('sq', num2str(ROI), '-onsetLatency(ms)AVG')};
% %     squareVariables(startingIndex+3) = {strcat('sq', num2str(ROI), '-latencyToPeak(ms)AVG')};
% %     squareVariables(startingIndex+4) = {strcat('sq', num2str(ROI), '-riseTime(ms)AVG')};
% %     squareVariables(startingIndex+5) = {strcat('sq', num2str(ROI), '-lightEvokedPeak(pA)STD')};
% %     squareVariables(startingIndex+6) = {strcat('sq', num2str(ROI), '-onsetLatency(ms)STD')};
% %     squareVariables(startingIndex+7) = {strcat('sq', num2str(ROI), '-latencyToPeak(ms)STD')};
% %     squareVariables(startingIndex+8) = {strcat('sq', num2str(ROI), '-riseTime(ms)STD')};
% % end
% % 
% % % label & save stored data
% % labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
% %     {'mouse', ...
% %     'date', ...
% %     'firstSweep', ...
% %     'lastSweep', ...
% %     'nSweeps', ...
% %     'lightChannel', ...
% %     'singleLightPulse(1)orTrain(0)', ...
% %     'lightPulseDur(s)', ... 
% %     'lightStimFreq(Hz)', ...
% %     'lightDur(s)', ...   
% %     'gridColumns', ...
% %     'gridRows', ...
% %     'inward(-1)ORoutward(1)', ...
% %     'baselineDurationInSeconds', ...
% %     'lightPulseAnalysisWindowInSeconds', ...
% %     'thresholdInDataPts', ...
% %     'rsTestPulseOnsetTime', ...
% %     'AVGseriesResistance(Mohm)', ...
% %     'MINseriesResistance(Mohm)', ...
% %     'MAXseriesResistance(Mohm)', ...
% %     'AVGbaselineCurrent(pA)', ...
% %     'MINbaselineCurrent(pA)', ...
% %     'MAXbaselineCurrent(pA)', ...
% %     squareVariables{:}, ...
% %     'cellImageFileNameDIC', ...
% %     'cellImageFileNameAlexa'});
% % writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
% % disp('I saved the square_by_square cell xls file')