% Align all detected action potentials and calculate action potential width

function apWidth(obj)

%% USER INPUT

% Saving CSV file
savefileto = 'D:\CORONAVIRUS DATA\From MATLAB';

% Illustrating APs
prePeak = 0.005;            % in seconds 
postPeak = 0.005;           % in seconds
baselineDuration = 0.002;   % in seconds

% Finding APs
peaksOrValleys = 'v';
highpassThreshold = 100;
lowpassThreshold = 1500;
MinPeakHeight = 10;
MinPeakDistance = 0.005;
discardedSweeps = 1;

% Analyzing AP
ddyValleyThreshold = 30;
ddyPeakThreshold = 30;

%% MAIN CODE

% getting info from file
mouseNumber = getMouseNumber(obj);
experimentDate = getExperimentDate(obj);
samplingFrequency = obj.header.Acquisition.SampleRate;

% finding sweep numbers from file name
[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj); 

numel(fieldnames(obj.sweeps));
obj.header.NSweepsPerRun;

% checking for incomplete sweeps and not analyzing incomplete sweeps - to
% avoid this error: "Index exceeds the number of array elements (0)".      
if numel(fieldnames(obj.sweeps)) <= obj.header.NSweepsPerRun  
    lastSweepNumber = firstSweepNumber + numel(fieldnames(obj.sweeps)) - discardedSweeps;
    allSweeps = firstSweepNumber:lastSweepNumber;
end 

% matrices that will be filled
xSubset = [];
ySubset = [];
yFilteredSubset = [];
xSubsetAll = [];
ySubsetAll = [];
yFilteredSubsetAll = [];
data = [];
nAPtotal = 0;

% Look for action potentials in all sweeps
for sweepNumber = allSweeps

    % Collect raw data from every sweep, channel 1
    [x,y] = obj.xy(sweepNumber, 1);
    
    % Bandpass filter data
    yFiltered = bandpass(y,[highpassThreshold lowpassThreshold],samplingFrequency);

    % Find action potentials - timestamps are stored in variable "locs"
    if peaksOrValleys == 'peaks'
        [pks,locs,w,p] = findpeaks(yFiltered,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
    else
        [pks,locs,w,p] = findpeaks(-yFiltered,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
    end

    % Store total number of APs found
    nAP = size(locs,1);
    
    % Store data surrounding each AP
    % Exclude first and last AP to avoid error - MATLAB will complain if
    % the required time window around each AP is beyond window of data
    % available.
    for iAP = 2:(nAP-1)
                
        % convert time points to data points
        xmin = round(samplingFrequency * (locs(iAP) - prePeak));
        xmax = round(samplingFrequency * (locs(iAP) + postPeak));
        xbaseline = round(xmin + samplingFrequency*baselineDuration);
        
        % make data rows  
        ySubset = transpose(y(xmin:xmax));
        yFilteredSubset = transpose(yFiltered(xmin:xmax));
        
        % calculate mean y during user-defined baseline
        ySubsetBaseline = mean(y(xmin:xbaseline));      
        
        % calculate baseline subtracted y
        ySubsetBaselineSubtracted = ySubset - ySubsetBaseline;
        
        % Store all data
        ySubsetAll = [ySubsetAll; ySubsetBaselineSubtracted];
        yFilteredSubsetAll = [yFilteredSubsetAll; yFilteredSubset];
        
    end   
    nAPtotal = nAPtotal + nAP;
    
end

% Calculate average AP
avgAP = mean(ySubsetAll);
avgAPfiltered = mean(yFilteredSubsetAll);

% Create x axis for plotting
xSubset = 1000*linspace(0,(prePeak + postPeak), 101);

% Find AP peak and valley
avgAPpeakInDataPoints = find(avgAP==max(avgAP));
avgAPvalleyInDataPoints = find(avgAP==min(avgAP));
avgAPpeakInMilliSeconds = xSubset(avgAPpeakInDataPoints);
avgAPvalleyInMilliSeconds = xSubset(avgAPvalleyInDataPoints);

% Find AP offset
% Assumes that AP valley precedes the AP peak
avgAPafterPeak = avgAP;
avgAPafterPeak(1:avgAPpeakInDataPoints) = [];
nExcludedDataPoints = size(avgAP, 2) - size(avgAPafterPeak, 2);
avgAPoffsetInDataPoints = find(round(avgAPafterPeak)==0, 1);
avgAPoffsetInMilliSeconds = xSubset(nExcludedDataPoints + avgAPoffsetInDataPoints);

% To avoid errors in the CSV in case matlab fails to find the offset:
if isempty(avgAPoffsetInMilliSeconds)
    avgAPoffsetInMilliSeconds = xSubset(end);
    avgAPoffsetInDataPoints = find(xSubset==max(xSubset)); 
    nExcludedDataPoints = 0;
end

% Find AP onset
% Assumes that AP valley precedes the AP peak
avgAPbeforeValley = avgAP;
avgAPbeforeValley(avgAPvalleyInDataPoints:end) = [];
avgAPonsetInDataPoints = find(round(avgAPbeforeValley)==0, 1, 'last');
avgAPonsetInMilliSeconds = xSubset(avgAPonsetInDataPoints);

% Calculating derivatives and creating xAxis for derivative plotting
xSubsetDropLast = xSubset;
xSubsetDropLast(end) = [];
dy = diff(avgAP)./diff(xSubset);
ddy = diff(dy)./diff(xSubsetDropLast);
xSubsetDropLastAgain = xSubsetDropLast;
xSubsetDropLastAgain(end) = [];

% Find AP onset based on 2nd derivative - first valley
[pks1,locs1,w,p] = findpeaks(-ddy,xSubsetDropLastAgain,'MinPeakHeight',ddyValleyThreshold);
ddyBasedOnset = locs1(1);         % in milli seconds
ddyFirstValley = pks1(1);
ddyBasedOnsetPt = round(ddyBasedOnset*(samplingFrequency/1000));

% Find AP offset based on 2nd derivative - last peak
[pks2,locs2,w,p] = findpeaks(ddy,xSubsetDropLastAgain,'MinPeakHeight',ddyPeakThreshold);
ddyBasedOffset = locs2(end);      % in milli seconds
ddyLastPeak = pks2(end);
ddyBasedOffsetPt = round(ddyBasedOffset*(samplingFrequency/1000));
ddyAfterLastPeak = ddy;
ddyAfterLastPeak(1:ddyBasedOffsetPt) = [];
ddyBasedOffset2Pt = ddyBasedOffsetPt + find(round(ddyAfterLastPeak)==0, 1)
ddyBasedOffset2 = xSubsetDropLastAgain(ddyBasedOffset2Pt)

% % Find AP onset based on 2nd derivative
% % Assumes that AP valley precedes the AP peak
% ddyPeakPt = find(ddy==max(ddy),1);
% ddyPrePeak = ddy;
% ddyPrePeak(ddyPeakPt:end) = [];
% ddyFirstValleyPt = find(ddy==min(ddyPrePeak),1); 
% 
% % Option 1 - onset is at ddy=0
% ddyPreFirstValley = ddy;
% ddyPreFirstValley(ddyFirstValleyPt:end) = [];
% ddyBasedOnsetPt = find(round(ddyPreFirstValley)==0, 1, 'last');
% ddyBasedOnsetMilliSec = xSubset(ddyBasedOnsetPt);
% 
% % Option 2 - onset is at a ddy valley
% % ddyBasedOnsetPt = ddyFirstValleyPt;
% % ddyBasedOnsetMilliSec = xSubset(ddyBasedOnsetPt);

% Calculate AP width
halfWidth = avgAPpeakInMilliSeconds - avgAPvalleyInMilliSeconds;
biphasicDuration = avgAPpeakInMilliSeconds - avgAPonsetInMilliSeconds;
duration = avgAPoffsetInMilliSeconds - avgAPonsetInMilliSeconds;
ddyBasedBiphasicDuration = avgAPpeakInMilliSeconds - ddyBasedOnset; 
ddyBasedDuration1 = avgAPoffsetInMilliSeconds - ddyBasedOnset;
ddyBasedDuration2 = ddyBasedOffset - ddyBasedOnset;


% Plot all APs and avg AP (raw, baseline subtracted)
figure('name', strcat(obj.file,' all APs Raw Baseline Subtracted'));     
hold on;
    plot(xSubset, ySubsetAll,'Color', [0.75, 0.75, 0.75, 0.5], 'LineWidth', 0.2);
    plot(xSubset, avgAP,'Color','black','LineWidth',1.5); 
    plot(xSubset(avgAPpeakInDataPoints), max(avgAP),'o','color', 'r');
    plot(xSubset(avgAPvalleyInDataPoints), min(avgAP),'o','color', 'b');
    plot(avgAPoffsetInMilliSeconds, avgAP(nExcludedDataPoints + avgAPoffsetInDataPoints),'o','color', 'g');
    plot(ddyBasedOffset, avgAP(ddyBasedOffsetPt), '*', 'color', 'g');
    plot(ddyBasedOffset2, avgAP(ddyBasedOffset2Pt), '+', 'color', 'g');
    plot(avgAPonsetInMilliSeconds, avgAP(avgAPonsetInDataPoints), 'o', 'color', 'y');
    plot(ddyBasedOnset, avgAP(ddyBasedOnsetPt), '*', 'color', 'y');
    xlabel('Time (ms)');
    ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
    title([obj.file ' all APs Raw Baseline Subtracted'],'Interpreter','none');
hold off;
movegui('east');

% Plot all APs and avg AP (bandpass filtered)
figure('name', strcat(obj.file,' all APs filtered'));     
hold on;
    plot(xSubset, yFilteredSubsetAll,'Color', [0.75, 0.75, 0.75, 0.5], 'LineWidth', 0.2);
    plot(xSubset, avgAPfiltered,'Color','black','LineWidth',1.5); 
    xlabel('Time (ms)');
    ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
    title([obj.file ' all APs filtered'],'Interpreter','none');
hold off;

% Plot the first (blue) and second (red) derivative of the avg
figure('name', strcat(obj.file,' avg Raw Baseline Subtracted Dvs'));
hold on;
    plot(xSubset, avgAP,'Color','black','LineWidth',1.5);
    plot(xSubsetDropLast, dy,'Color', 'b', 'LineWidth', 1);
    plot(xSubsetDropLastAgain, ddy,'Color', 'r', 'LineWidth', 1);
    plot(ddyBasedOffset, ddy(ddyBasedOffsetPt), '*', 'color', 'g');
    plot(ddyBasedOffset2, ddy(ddyBasedOffset2Pt), '+', 'color', 'g');
    plot(ddyBasedOnset, ddy(ddyBasedOnsetPt), '*', 'color', 'y');
    xlabel('Time (ms)');
    ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
    title([obj.file ' avg Raw Baseline Subtracted Dvs'],'Interpreter','none');
hold off;
movegui('west');

% Store all data
data = [mouseNumber, ...
    experimentDate, ...
    firstSweepNumber, ...
    lastSweepNumber, ...
    prePeak, ...
    postPeak, ...
    baselineDuration, ...
    highpassThreshold, ...
    lowpassThreshold, ...
    MinPeakHeight, ...
    MinPeakDistance, ...
    discardedSweeps, ...
    nAPtotal, ...
    avgAPonsetInMilliSeconds, ...
    ddyBasedOnset, ...
    avgAPvalleyInMilliSeconds, ...
    avgAPpeakInMilliSeconds, ...
    avgAPoffsetInMilliSeconds, ...
    ddyBasedOffset, ...
    halfWidth, ...
    biphasicDuration, ...
    duration, ...
    ddyBasedBiphasicDuration, ...
    ddyBasedDuration1, ...
    ddyBasedDuration2];    

% save csv file with data 
filename = strcat(obj.file, " - AP width");
fulldirectory = strcat(savefileto,'\',filename,'.csv');        
dataInCellFormat = {};
dataInCellFormat = num2cell(data);
labeledData = cell2table(dataInCellFormat,'VariableNames',...
    {'mouse', 'date', 'firstSweep', 'lastSweep', ...
    'prePeak', ...
    'postPeak', ...
    'baselineDuration', ...
    'highpassThreshold', ...
    'lowpassThreshold', ...
    'MinPeakHeight', ...
    'MinPeakDistance', ...
    'discardedSweeps', ...
    'nAPtotal', ...
    'avgAPonsetInMilliSeconds', ...
    'ddyBasedOnsetMilliSec', ...
    'avgAPvalleyInMilliSeconds', ...
    'avgAPpeakInMilliSeconds', ...
    'avgAPoffsetInMilliSeconds', ...
    'ddyBasedOffsetMilliSec', ...
    'halfWidth', ...
    'biphasicDuration', ...
    'duration', ...
    'ddyBasedBiphasicDuration', ...
    'ddyBasedDuration1', ...
    'ddyBasedDuration2'})
writetable(labeledData,fulldirectory);
disp('I saved the CSV, you save the pics')
   
end    