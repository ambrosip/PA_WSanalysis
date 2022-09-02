%{ 
DOCUMENTATION
Created: 2022 08 30
Edited Last: 2022 09 01
Works? Yes
Author: PA

This function is used to analyze light-evoked PSCs in sCRACM experiments using the polygon.
This function was derived from firing_vs_light_polygon, but went through A LOT of rewritting.

INPUTS explained:
    - gridColumns: number of columns in the polygon ROI grid

    - gridRows: number of rows in the polygon ROI grid

    - discardedSweeps: specific sweeps that should NOT be analyzed due to
    artifacts. If I want to discard sweep 0024, just write 24.

    - discardedSweepsFromEnd: number of sweeps from the end that will not be
    analyzed. Used to remove incomplete sweeps or sweeps with artifacts. 

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
    acommodate the identification of smaller PSCs)

    - rsTestPulseOnsetTime: onset of voltage step used to calculate series
    resistance.

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
    rsTestPulseOnsetTime = 1;

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

function psc_vs_light_polygon(obj)
%%  USER INPUT ============================================================

% Affects data analysis - Organizing data by o-stim grid
gridColumns = 5;
gridRows = 5;

% Affects data analysis - Finding/quantifyting oIPSCs
discardedSweeps = [];
discardedSweepsFromEnd = 0;
lightChannel = 4;
singleLightPulse = 1; 
inwardORoutward = -1;               % 1 (positive) is outward; -1 (negative) in inward
baselineDurationInSeconds = 0.01;
lightPulseAnalysisWindowInSeconds = 0.02;
thresholdInDataPts = 5;             % ALERT! Changed from 10 to 5
rsTestPulseOnsetTime = 1;

% Affects data display: 
ymin = -1500;           %-2050
ymax = 600;             %50
cellImageFileNameDIC = 's3c1_z1_dic.tif';
cellImageFileNameAlexa = 's3c1_z1_647.tif';
cellImageDir = 'D:\NU server\Priscilla - BACKUP 20200319\Ephys\2022\20220727 m075 dls loop';

% Affects data saving:
savefileto = 'D:\Temp\From MATLAB 2022 09 31 psc better';

% % Affects oIPSC decay fit - Not currently implemented in polygon code
% bGuess = 1;             % -1000     % 1
% fastTauGuess = 0.01;    % in seconds
% slowTauGuess = 0.1;     % in seconds


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

% calculating variables based on user input
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
lightEvokedCurrentsAmp = [];
allLightEvokedResponseOnsetLatencyInMilliSeconds = [];
allRs = [];
yBaselineSubAll = [];
lightEvokedCurrentsAllSweeps = [];
allLightEvokedResponseOnsetLatencyInMilliSecondsAllSweeps = [];
allLightEvokedResponsePeakLatencyInMilliSecondsAllSweeps = [];
allTimeTo10percentOfPeakInMilliSecondsAllSweeps = [];
allTimeTo90percentOfPeakInMilliSecondsAllSweeps = [];
allRiseTimeInMilliSecondsAllSweeps = [];
baselineCurrentAll = [];
data = [];
isDiscardedBySweep = [];


%% SWEEP BY SWEEP ANALYSIS ===================================================================================

% get data from all sweeps in file
for sweepNumber = allSweeps
    
    % get light stim data
    [xch2,ych2] = obj.xy(sweepNumber, lightChannel); 

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
%         xmax = lightOnsetTime + lightDur + 2*baselineDurationInSeconds; 
%         xmax = lightOnsetTime + lightDur + baselineDurationInSeconds; 
%         xmax = lightOnsetTime+0.2;        
    end    
    
    % set xmin for plotting later
    xmin = lightOnsetTime - baselineDurationInSeconds;
    %----------------------------------------------------------------------
    
    
    % ignore discarded sweeps, but get raw data from other sweeps
    if ismember(sweepNumber, discardedSweeps)        
        isDiscarded = 1;
        x = NaN(sweepDurationInDataPts, 1);
        y = NaN(sweepDurationInDataPts, 1);
    else         
        isDiscarded = 0;
        [x,y] = obj.xy(sweepNumber, 1);  
    end
    
    % store info about which sweeps were discarded
    isDiscardedBySweep = [isDiscardedBySweep; isDiscarded];
    %----------------------------------------------------------------------
    
    
    % baseline subtraction
    baselineStart = lightPulseStart(1) - baselineDurationInDataPts;
    yBaselineSub = y-mean(y(baselineStart:lightPulseStart(1)));
    baselineCurrent = mean(y(baselineStart:lightPulseStart(1)));
    
    % saving data for niceplot
    % y data for each sweep is in a column
    yBaselineSubAll = [yBaselineSubAll, yBaselineSub]; 
    baselineCurrentAll = [baselineCurrentAll, baselineCurrent];
    %----------------------------------------------------------------------
    
    
    % calculating series resistance
    rsBaselineCurrent = mean(y(rsBaselineDataPointInterval));
    rsTransientCurrent = min(y(rsFirstTransientDataPointInterval));
    dCurrent = rsTransientCurrent-rsBaselineCurrent;
    dVoltage = -5;
    seriesResistance = 1000*dVoltage/dCurrent; %mV/pA equals Gohm
    
    % put series resistance values from each sweep into a different column
    allRs = [allRs, seriesResistance];
    %----------------------------------------------------------------------
    
    
    % clean up matrices that will be used in the next loop
    lightEvokedCurrentsAmp = [];
    lightEvokedCurrentsLoc = [];
    allLightEvokedResponseOnsetLatencyInMilliSeconds = [];
    allLightEvokedResponsePeakLatencyInMilliSeconds = [];
    allTimeTo10percentOfPeakInMilliSeconds = [];
    allTimeTo90percentOfPeakInMilliSeconds = [];
    allRiseTimeInMilliSeconds = [];
    
    % loop through all light pulses in the train
    % this code is overkill here cuz there is only one light pulse but whatever
    for pulseOnset = lightPulseStart.'
        
        afterLightDataPoint = pulseOnset + lightPulseAnalysisWindowInDataPts;
            
        % get amplitude and location (index) of lightEvokedCurrents
        % outward current is +1
        % inward current is 0 or -1
        if inwardORoutward == 1 
            [lightEvokedCurrentAmp, lightEvokedCurrentLoc] = max(yBaselineSub(pulseOnset:afterLightDataPoint));
            timeTo10percentOfPeakInDataPoints = find(yBaselineSub(pulseOnset:afterLightDataPoint) >= 0.1 * lightEvokedCurrentAmp, 1);
            timeTo90percentOfPeakInDataPoints = find(yBaselineSub(pulseOnset:afterLightDataPoint) >= 0.9 * lightEvokedCurrentAmp, 1);
        else
            [lightEvokedCurrentAmp, lightEvokedCurrentLoc] = min(yBaselineSub(pulseOnset:afterLightDataPoint));
            timeTo10percentOfPeakInDataPoints = find(yBaselineSub(pulseOnset:afterLightDataPoint) <= 0.1 * lightEvokedCurrentAmp, 1);
            timeTo90percentOfPeakInDataPoints = find(yBaselineSub(pulseOnset:afterLightDataPoint) <= 0.9 * lightEvokedCurrentAmp, 1);
        end
        
        % put all currents from a particular sweep in a row (each column is
        % a pulse)
        lightEvokedCurrentsAmp = [lightEvokedCurrentsAmp, lightEvokedCurrentAmp];
        lightEvokedCurrentsLoc = [lightEvokedCurrentsLoc, lightEvokedCurrentLoc];
               
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
        
        % put all latencies from a particular sweep in a row
        allLightEvokedResponseOnsetLatencyInMilliSeconds = [allLightEvokedResponseOnsetLatencyInMilliSeconds, lightEvokedResponseOnsetLatencyInMilliSeconds];        
        allLightEvokedResponsePeakLatencyInMilliSeconds = [allLightEvokedResponsePeakLatencyInMilliSeconds, lightEvokedResponsePeakLatencyInMilliSeconds];        
        allTimeTo10percentOfPeakInMilliSeconds = [allTimeTo10percentOfPeakInMilliSeconds, timeTo10percentOfPeakInMilliSeconds];
        allTimeTo90percentOfPeakInMilliSeconds = [allTimeTo90percentOfPeakInMilliSeconds, timeTo90percentOfPeakInMilliSeconds];
        allRiseTimeInMilliSeconds = [allRiseTimeInMilliSeconds, riseTimeInMilliSeconds];
           
    end
    
    % store all currents from a cell (each column is a pulse; each row is a sweep)
    lightEvokedCurrentsAllSweeps = [lightEvokedCurrentsAllSweeps; lightEvokedCurrentsAmp];
    allLightEvokedResponseOnsetLatencyInMilliSecondsAllSweeps = [allLightEvokedResponseOnsetLatencyInMilliSecondsAllSweeps; allLightEvokedResponseOnsetLatencyInMilliSeconds];
    allLightEvokedResponsePeakLatencyInMilliSecondsAllSweeps = [allLightEvokedResponsePeakLatencyInMilliSecondsAllSweeps; allLightEvokedResponsePeakLatencyInMilliSeconds];
    allTimeTo10percentOfPeakInMilliSecondsAllSweeps = [allTimeTo10percentOfPeakInMilliSecondsAllSweeps; allTimeTo10percentOfPeakInMilliSeconds];
    allTimeTo90percentOfPeakInMilliSecondsAllSweeps = [allTimeTo90percentOfPeakInMilliSecondsAllSweeps; allTimeTo90percentOfPeakInMilliSeconds];
    allRiseTimeInMilliSecondsAllSweeps = [allRiseTimeInMilliSecondsAllSweeps; allRiseTimeInMilliSeconds];    
    
    % store key data that will be exported    
    % each sweep is a row
    data = [data; ...
        mouseNumber, ...
        experimentDate, ...
        sweepNumber, ...
        isDiscarded, ...
        discardedSweepsFromEnd, ...           
        lightChannel, ...
        singleLightPulse, ...
        stimDur, ... 
        stimFreq, ...
        lightDur, ...
        gridColumns, ...
        gridRows, ...
        sweepsPerSquare, ...
        inwardORoutward, ...
        baselineDurationInSeconds, ...
        lightPulseAnalysisWindowInSeconds, ...
        thresholdInDataPts, ...
        rsTestPulseOnsetTime, ...
        seriesResistance, ...
        baselineCurrent];             
end

% store key data that will be exported
data = [data, ... 
    lightEvokedCurrentsAllSweeps, ...
    allLightEvokedResponseOnsetLatencyInMilliSecondsAllSweeps, ...
    allLightEvokedResponsePeakLatencyInMilliSecondsAllSweeps, ...
    allTimeTo10percentOfPeakInMilliSecondsAllSweeps, ...
    allTimeTo90percentOfPeakInMilliSecondsAllSweeps, ...
    allRiseTimeInMilliSecondsAllSweeps];


%% CELL ANALYSIS - POLYGON SPECIFIC CODE - data for each square ==============================================

% Store cell-specific data
dataCell = [mouseNumber, ...
    experimentDate, ...
    firstSweepNumber, ...
    lastSweepNumber, ...
    length(allSweeps) - length(discardedSweeps), ...     
    lightChannel, ...
    singleLightPulse, ...
    stimDur, ... 
    stimFreq, ...
    lightDur, ...
    gridColumns, ...
    gridRows, ...
    inwardORoutward, ...
    baselineDurationInSeconds, ...
    lightPulseAnalysisWindowInSeconds, ...
    thresholdInDataPts, ...
    rsTestPulseOnsetTime, ...
    mean(allRs), ...
    min(allRs), ...
    max(allRs), ...
    mean(baselineCurrentAll), ...
    min(baselineCurrentAll), ...
    max(baselineCurrentAll)];

% Reintroduce these bois up there if you ever want to analyze the decay.
% And don't forget to reintroduce them in the "EXPORTING" segment as well as
% variable names.
%     bGuess, ...
%     fastTauGuess, ...
%     slowTauGuess, ...

% Assigning sweeps to each o-stim area (square or ROI)
% Each row in sweepIDperSquare correspond to a square
% Each column in sweepIDperSquare corresponds to a repeat sweep in that
% square.
totalSquares = gridColumns * gridRows;
relativeSweepNumberPerSquare = []; % from 1 to totalSquares
absoluteSweepNumberPerSquare = []; % from 1st sweep to last sweep number
dataSquare = [];
yBaselineSubAllMeanBySquare = [];
dataCellSingleRow = dataCell;

% keep track of sweeps per square taking into consideration any
% discarded sweeps. 
sweepsPerSquarePerSquareAll = [];

% each row is a square
for row=[1:totalSquares]
    
    % create/erase arrays that will be filled for each square
    dataSubsetForMeanAndSD = [];
    dataSubsetForNicePlot = [];
    
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
            % yBaselineSubAll is organized with 1 sweep per column
            % to get the full data from one sweep, I need all elements in a
            % particular column (:,column#)
            % dataSubsetForNicePlot is organized with 1 sweep per column as
            % well (comma separates new entries)
            dataSubsetForNicePlot = [dataSubsetForNicePlot, ...
                                        yBaselineSubAll(:,relativeSweepNumber)];
            
            % ALERT need to check the organization of these matrixes
            % data from each sweep is in a row
            dataSubsetForMeanAndSD = [dataSubsetForMeanAndSD; ...
                                        lightEvokedCurrentsAllSweeps(relativeSweepNumber), ...
                                        allLightEvokedResponseOnsetLatencyInMilliSecondsAllSweeps(relativeSweepNumber), ...
                                        allLightEvokedResponsePeakLatencyInMilliSecondsAllSweeps(relativeSweepNumber), ...
                                        allRiseTimeInMilliSecondsAllSweeps(relativeSweepNumber)];                                                        
        end         
    end
       
    % store number of sweeps analyzed for each square
    sweepsPerSquarePerSquare = size(dataSubsetForMeanAndSD,1);
    sweepsPerSquarePerSquareAll = [sweepsPerSquarePerSquareAll; sweepsPerSquarePerSquare];
    
    % calculate mean across sweeps for each square
    yBaselineSubAllMean = mean(dataSubsetForNicePlot,2,'omitnan');
    
    % store mean by square (data from each square in a different column)
    yBaselineSubAllMeanBySquare = [yBaselineSubAllMeanBySquare, yBaselineSubAllMean];    
        
    % store data for each square
    % each square is a row
    dataSquare = [dataSquare;...
        dataCell,...                             % all the cell info
        row,...                                  % square number
        absoluteSweepNumberPerSquare(row,:),...  % all sweeps in this square        
        sweepsPerSquarePerSquare,...
        mean(dataSubsetForMeanAndSD, 'omitnan'),...
        std(dataSubsetForMeanAndSD, 'omitnan')];
                
    % store square-by-square data together with cell data
    dataCellSingleRow = [dataCellSingleRow, ...
                    sweepsPerSquarePerSquare, ...
                    mean(dataSubsetForMeanAndSD, 'omitnan'), ...
                    std(dataSubsetForMeanAndSD, 'omitnan')];                 
                
    % alternative code to add dataCell to dataSquare
% % concatenate dataSquare with dataCell
% % use repmat to copy dataCell into multiple rows - as many rows as there
% % are squares
% dataCellMultipleRows = [repmat(dataCell,totalSquares,1), dataSquare];
    
    % commented out all average signal analysis for the polygon dataset on
    % 2022 09 01 because it is not reasonable to analyze the average of 2-3
    % sweeps per square. It is misleading.
% %     % analyze the mean signal per square
% %     % get amplitude and location (index) of lightEvokedCurrents
% %     % for rise time calculation, find index of 10% and 90% of peak
% %     % look for a threshold cross instead of an exact match - find(a<=b) instead of find(a==b) 
% %     % outward current is +1
% %     % inward current is 0 or -1
% %     % ASSUMPTION ALERT: Assuming that all sweeps have the same, single o-stim
% %     pulseOnset = lightPulseStart(1);
% %     afterLightDataPoint = pulseOnset + lightPulseAnalysisWindowInDataPts;
% %     if inwardORoutward == 1 
% %         [lightEvokedCurrentAmp, lightEvokedCurrentLoc] = max(yBaselineSubAllMean(pulseOnset:afterLightDataPoint));
% %         latencyTo10percentOfPeakInDataPoints = find(yBaselineSubAllMean(pulseOnset:afterLightDataPoint) >= 0.1 * lightEvokedCurrentAmp, 1) + pulseOnset;
% %         latencyTo90percentOfPeakInDataPoints = find(yBaselineSubAllMean(pulseOnset:afterLightDataPoint) >= 0.9 * lightEvokedCurrentAmp, 1) + pulseOnset;
% %     else
% %         [lightEvokedCurrentAmp, lightEvokedCurrentLoc] = min(yBaselineSubAllMean(pulseOnset:afterLightDataPoint));
% %         latencyTo10percentOfPeakInDataPoints = find(yBaselineSubAllMean(pulseOnset:afterLightDataPoint) <= 0.1 * lightEvokedCurrentAmp, 1) + pulseOnset;
% %         latencyTo90percentOfPeakInDataPoints = find(yBaselineSubAllMean(pulseOnset:afterLightDataPoint) <= 0.9 * lightEvokedCurrentAmp, 1) + pulseOnset;
% %     end    
% % 
% %     % adjust latency to peak Loc
% %     lightEvokedCurrentLoc = lightEvokedCurrentLoc + pulseOnset;
% % 
% %     % find onset latency
% %     onsetLatencyInDataPoints = min(find(conv2(sign(diff(yBaselineSubAllMean(pulseOnset:afterLightDataPoint))), zeros(thresholdInDataPts,1)+(inwardORoutward), 'valid')==thresholdInDataPts));
% %     onsetLatencyInDataPoints = onsetLatencyInDataPoints + pulseOnset;
    
    % decay analysis - commented out on 2022 09 01 because the polygon
    % dataset is not the best for this kind of calculation - many squares
    % do not have oIPSCs, so decay calculation is meaningless (AND SLOW)
% %     % for decay tau calculation, find index of return to baseline after peak
% %     latencyToZeroAfterPeakInDataPoints = find(round(yBaselineSubAllMean(lightEvokedCurrentLoc:end)) >= 0, 1) + lightEvokedCurrentLoc
% % 
% %     % preparing to fit double-term exponential model to the decay data - create y subset
% %     xminFit = lightEvokedCurrentLoc;
% %     xmaxFit = latencyToZeroAfterPeakInDataPoints;
% %     yBaselineSubAllMeanSubset = yBaselineSubAllMean;
% %     yBaselineSubAllMeanSubset(latencyToZeroAfterPeakInDataPoints:end) = [];
% %     yBaselineSubAllMeanSubset(1:lightEvokedCurrentLoc) = [];
% % 
% %     % preparing to fit double-term exponential model to the decay data - create x subset
% %     xSubset = x;
% %     xSubset(latencyToZeroAfterPeakInDataPoints:end) = [];
% %     xSubset(1:lightEvokedCurrentLoc) = [];
% %     xSubset = round(xSubset,6) - round(lightEvokedCurrentLoc/samplingFrequency,6);
% % 
% %     % selecting a double exponential for decay fitting
% %     g = fittype('a*exp(-x/fastTau) + b*exp(-x/slowTau)');
% % 
% %     % fit double-term exponential model to the decay data 
% %     % If you don't give MATLAB a StartPoint, it fails horribly
% %     % StartPoint order: [a, b, fastTau, slowTau]
% %     % You can check the order by using this line:
% %     % coefficientNames = coeffnames(g)
% %     % I guessed StartPoint values by manually adjusting a double exponential
% %     % The single exponential was NOT giving me a good fit
% %     [decayFit,gof,output] = fit(xSubset,yBaselineSubAllMeanSubset,g,'StartPoint',[lightEvokedCurrentAmp, bGuess, fastTauGuess, slowTauGuess]);
% % 
% %     % save coefficients a and b
% %     decayFitA = decayFit.a;
% %     decayFitB = decayFit.b;
% % 
% %     % convert s to ms
% %     decayFitFastTau = decayFit.fastTau * 1000
% %     decayFitSlowTau = decayFit.slowTau * 1000     
% % 
% %     % % For troubleshooting and estimating StartPoints:
% %     % figure;
% %     % plot(decayFit);
% %     % figure;
% %     % hold on;
% %     % plot(xSubset+round(lightEvokedCurrentLoc/samplingFrequency,6), decayFit(xSubset))
% %     % plot(x, yBaselineSubAllMean);
% %     % plot(xSubset+round(lightEvokedCurrentLoc/samplingFrequency,6), lightEvokedCurrentAmp*exp(-xSubset/0.05));
% %     % hold off;

    % continue to comment out average analysis
% %     % convert data points to ms
% %     onsetLatencyInMilliSeconds = 1000 * (onsetLatencyInDataPoints - pulseOnset) / samplingFrequency;
% %     latencyToPeakInMilliSeconds = 1000 * (lightEvokedCurrentLoc - pulseOnset) / samplingFrequency;
% %     latencyToZeroAfterPeakInMilliSeconds = 1000 * (latencyToZeroAfterPeakInDataPoints - pulseOnset) / samplingFrequency;
% %     latencyTo10percentOfPeakInMilliSeconds = 1000 * (latencyTo10percentOfPeakInDataPoints - pulseOnset) / samplingFrequency;
% %     latencyTo90percentOfPeakInMilliSeconds = 1000 * (latencyTo90percentOfPeakInDataPoints - pulseOnset) / samplingFrequency;
% % 
% %     % check if onset latency could be calculated to avoid errors
% %     if isempty(onsetLatencyInDataPoints)
% %         onsetLatencyInMilliSeconds = NaN;
% %         onsetLatencyInDataPoints = NaN;
% %     end
% % 
% %     % calculate riseTime in MilliSeconds
% %     riseTimeInMilliSeconds = latencyTo90percentOfPeakInMilliSeconds - latencyTo10percentOfPeakInMilliSeconds;
  
end


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
figure('name', strcat(fileName, '_', analysisDate, ' - psc_vs_light_polygon - DIC image'));
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


%% PLOT 2 - cropped DIC cell image with grid ===============================================================================
% make sure you're getting the image taken with zoom = 1
% concatenate strings from user input to get full path to figure file
cellImageFileDir = strcat(cellImageDir,'\',cellImageFileNameDIC);

% import image
cellImage = imread(cellImageFileDir);

% crop image according to the polygon's mirror calibration
croppedImage = imcrop(cellImage, [1,100,1376,873]);
figure('name', strcat(fileName, '_', analysisDate, ' - psc_vs_light_polygon - DIC image grid'));
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
set(gcf,'InnerPosition',[0 maxHeight-innerHeight innerWidth innerHeight]);


%% PLOT 3 - cropped Alexa cell image with grid ===============================================================================
% make sure you're getting the image taken with zoom = 1
% concatenate strings from user input to get full path to figure file
cellImageFileDir = strcat(cellImageDir,'\',cellImageFileNameAlexa);

% import image
cellImage = imread(cellImageFileDir);

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


%% PLOT 4 - cropped Alexa cell image without grid ===============================================================================
% make sure you're getting the image taken with zoom = 1
% concatenate strings from user input to get full path to figure file
cellImageFileDir = strcat(cellImageDir,'\',cellImageFileNameAlexa);

% import image
cellImage = imread(cellImageFileDir);

% invert alexa image so that black = white
% I wanted to do this anyway, but MATLAB also forced my hand. When I save
% images with my saveAllFigs function, MATLAB turns all white annotations
% (text and lines) from white to black, rendering my scale bar useless.
invertedImage = imcomplement(cellImage);

% crop image according to the polygon's mirror calibration
croppedImage = imcrop(invertedImage, [1,100,1376,873]);
figure('name', strcat(fileName, '_', analysisDate, ' - psc_vs_light_polygon - Alexa image'));
hold on;
imshow(croppedImage, 'Border', 'tight');

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
plot(allSweeps, allRs,'-o');
% plot lines marking 30% increase and 30% decrese in Rs compared to first
% test pulse
line([allSweeps(1) allSweeps(end)],[allRs(1)*0.7, allRs(1)*0.7],'Color','black','LineStyle','--')
line([allSweeps(1) allSweeps(end)],[allRs(1)*1.3, allRs(1)*1.3],'Color','black','LineStyle','--')
axis([allSweeps(1) inf 0 60])
ylabel('Rs (M\Omega)');
xlabel('Sweeps');
title([obj.file ' rs'],'Interpreter','none');
movegui('southeast');


%% PLOT 6 - Tiled oIPSC amplitude =========================================================================

% create figure & name it
figure('name', strcat(fileName, '_', analysisDate, ' - psc_vs_light_polygon - tiled niceplots'));
t = tiledlayout(gridRows, gridColumns);

% plotting histogram plot 
for square = [1:totalSquares]
    nexttile
    hold on;
    
    % plot individual sweeps
    for sweep = [relativeSweepNumberPerSquare(square,:)]
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
        line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-inwardORoutward*(ymax),-inwardORoutward*(ymax)],'Color',[0 0.4470 0.7410],'LineWidth',10)
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
    if square == gridRows * gridColumns
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
plot(allSweeps, baselineCurrentAll,'-o');
axis([allSweeps(1) inf -500 500])
ylabel('Subtracted Baseline Current (pA)');
xlabel('Sweeps');
title([obj.file ' baseline current'],'Interpreter','none');
movegui('southwest');


%% EXPORTING XLS files ==========================================

% create cell array with strings for naming the amplitude and latency of the oPSC for each light pulse
oPSCvariableNamesAmplitude = cell(1,length(lightEvokedCurrentsAmp));
oPSCvariableNamesOnsetLatency = cell(1,length(allLightEvokedResponseOnsetLatencyInMilliSeconds));
oPSCvariableNamesPeakLatency = cell(1,length(allLightEvokedResponsePeakLatencyInMilliSeconds));
oPSCvariableNames10percentLatency = cell(1,length(allTimeTo10percentOfPeakInMilliSeconds));
oPSCvariableNames90percentLatency = cell(1,length(allTimeTo90percentOfPeakInMilliSeconds));
oPSCvariableNamesRiseTime = cell(1,length(allRiseTimeInMilliSeconds));
for lightPulse = 1:length(lightEvokedCurrentsAmp)
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
    'isDiscarded(1)orNot(0)', ...
    'discardedSweepsFromEnd', ...
    'lightChannel', ...
    'singleLightPulse(1)orTrain(0)', ...
    'lightPulseDur(s)', ... 
    'lightStimFreq(Hz)', ...
    'lightDur(s)', ...   
    'gridColumns', ...
    'gridRows', ...
    'sweepsPerSquare', ...    
    'inward(-1)ORoutward(1)', ...
    'baselineDurationInSeconds', ...
    'lightPulseAnalysisWindowInSeconds', ...
    'thresholdInDataPts', ...
    'rsTestPulseOnsetTime', ...
    'seriesResistance(Mohm)', ...
    'baselineCurrent(pA)', ...
    oPSCvariableNamesAmplitude{:}, ...
    oPSCvariableNamesOnsetLatency{:}, ...
    oPSCvariableNamesPeakLatency{:}, ...
    oPSCvariableNames10percentLatency{:}, ... 
    oPSCvariableNames90percentLatency{:}, ...
    oPSCvariableNamesRiseTime{:}});

writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the sweep_by_sweep xls file')
%--------------------------------------------------------------------------


% warn user about current code limitations
% cuz I don't want to deal with a gigantic csv containing multiple light
% pulses for multiple squares!
if singleLightPulse == 0
    disp('only analysing the first oIPSC in each sweep bro, look out if you got more')
end    

% store key square-by-square data - each row is a square
filename = strcat(fileName, '_', analysisDate, " - psc_vs_light_polygon - square_by_square");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(dataSquare);
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
    'lightChannel', ...
    'singleLightPulse(1)orTrain(0)', ...
    'lightPulseDur(s)', ... 
    'lightStimFreq(Hz)', ...
    'lightDur(s)', ...   
    'gridColumns', ...
    'gridRows', ...
    'inward(-1)ORoutward(1)', ...
    'baselineDurationInSeconds', ...
    'lightPulseAnalysisWindowInSeconds', ...
    'thresholdInDataPts', ...
    'rsTestPulseOnsetTime', ...
    'AVGseriesResistance(Mohm)', ...
    'MINseriesResistance(Mohm)', ...
    'MAXseriesResistance(Mohm)', ...
    'AVGbaselineCurrent(pA)', ...
    'MINbaselineCurrent(pA)', ...
    'MAXbaselineCurrent(pA)', ...
    'square', ...
    headings{:}, ...
    'sweepsPerSquare', ...
    'AVGlightEvokedPeakCurrentAmp(pA)', ...
    'AVGonsetLatency(ms)', ...
    'AVGlatencyToPeak(ms)', ...
    'AVGriseTime(ms)', ...  
    'STDlightEvokedPeakCurrentAmp(pA)', ...
    'STDonsetLatency(ms)', ...
    'STDlatencyToPeak(ms)', ...
    'STDriseTime(ms)'});
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the square_by_square xls file, master')
%--------------------------------------------------------------------------


% store cell data with square-by-square info in a single row
filename = strcat(fileName, '_', analysisDate, " - psc_vs_light_polygon - square_by_square cell");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(dataCellSingleRow);
dataInCellFormat = [dataInCellFormat, cellImageFileNameDIC, cellImageFileNameAlexa];

% create cell array with strings for naming square-by-square data
% calculate the number of unique variables per square
numberOfVariablesPerSquare = size([sweepsPerSquarePerSquare, ...
                            mean(dataSubsetForMeanAndSD, 'omitnan'), ...
                            std(dataSubsetForMeanAndSD, 'omitnan')], 2);
% make a cell array with 1 row and as many columns as the number of unique
% variables we have (totalSquares*numberOfVariablesPerSquare)
squareVariables = cell(1,totalSquares*numberOfVariablesPerSquare);
% name each variable with square number and correct label
for square = [1:totalSquares]
    startingIndex = square*numberOfVariablesPerSquare+1 - numberOfVariablesPerSquare;
    squareVariables(startingIndex) = {strcat('sq', num2str(square), '-sweepsPerSquare')};
    squareVariables(startingIndex+1) = {strcat('sq', num2str(square), '-lightEvokedPeak(pA)AVG')};
    squareVariables(startingIndex+2) = {strcat('sq', num2str(square), '-onsetLatency(ms)AVG')};
    squareVariables(startingIndex+3) = {strcat('sq', num2str(square), '-latencyToPeak(ms)AVG')};
    squareVariables(startingIndex+4) = {strcat('sq', num2str(square), '-riseTime(ms)AVG')};
    squareVariables(startingIndex+5) = {strcat('sq', num2str(square), '-lightEvokedPeak(pA)STD')};
    squareVariables(startingIndex+6) = {strcat('sq', num2str(square), '-onsetLatency(ms)STD')};
    squareVariables(startingIndex+7) = {strcat('sq', num2str(square), '-latencyToPeak(ms)STD')};
    squareVariables(startingIndex+8) = {strcat('sq', num2str(square), '-riseTime(ms)STD')};
end

% label & save stored data
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
    {'mouse', ...
    'date', ...
    'firstSweep', ...
    'lastSweep', ...
    'nSweeps', ...
    'lightChannel', ...
    'singleLightPulse(1)orTrain(0)', ...
    'lightPulseDur(s)', ... 
    'lightStimFreq(Hz)', ...
    'lightDur(s)', ...   
    'gridColumns', ...
    'gridRows', ...
    'inward(-1)ORoutward(1)', ...
    'baselineDurationInSeconds', ...
    'lightPulseAnalysisWindowInSeconds', ...
    'thresholdInDataPts', ...
    'rsTestPulseOnsetTime', ...
    'AVGseriesResistance(Mohm)', ...
    'MINseriesResistance(Mohm)', ...
    'MAXseriesResistance(Mohm)', ...
    'AVGbaselineCurrent(pA)', ...
    'MINbaselineCurrent(pA)', ...
    'MAXbaselineCurrent(pA)', ...
    squareVariables{:}, ...
    'cellImageFileNameDIC', ...
    'cellImageFileNameAlexa'});
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the square_by_square cell xls file')

end