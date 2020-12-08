% Align all detected action potentials and calculate action potential width

function apWidthDdy(obj)

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
MinPeakHeight = 15;
MinPeakDistance = 0.05;
discardedSweeps = 5;

% Analyzing AP
ddyValleyThreshold = 30;
ddyPeakThreshold = 30;

%% Code warm-up
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

%% Look for action potentials in all sweeps
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

%% Find AP peak and valley
avgAPpeakInDataPoints = find(avgAP==max(avgAP));
avgAPvalleyInDataPoints = find(avgAP==min(avgAP));
avgAPpeakInMilliSeconds = xSubset(avgAPpeakInDataPoints);
avgAPvalleyInMilliSeconds = xSubset(avgAPvalleyInDataPoints);


%% Calculating derivatives and creating xAxis to plot derivatives 
% For each derivative, the xAxis length decreases by 1 point.
xForDy = xSubset;
xForDy(end) = [];   % remove last point from xAxis
dy = diff(avgAP)./diff(xSubset);
ddy = diff(dy)./diff(xForDy);
xForDdy = xForDy;
xForDdy(end) = [];  % remove last point from xAxis


%% Find AP onset based on 2nd derivative - ddy == min of first valley
[pks1,locs1,w,p] = findpeaks(-ddy,xForDdy,'MinPeakHeight',ddyValleyThreshold);
ddyBasedOnsetInMilliSeconds = locs1(1);       
ddyBasedOnsetInDataPoints = round(ddyBasedOnsetInMilliSeconds*(samplingFrequency/1000));


%% Find AP offset based on 2nd derivative - ddy == zero after last peak
[pks2,locs2,w,p] = findpeaks(ddy,xForDdy,'MinPeakHeight',ddyPeakThreshold);
ddyLastPeakInMilliSeconds = locs2(end);     
ddyLastPeakInDataPoints = round(ddyLastPeakInMilliSeconds*(samplingFrequency/1000));
ddyLastValleyInMilliSeconds = locs1(end);     
ddyLastValleyInDataPoints = round(ddyLastValleyInMilliSeconds*(samplingFrequency/1000));

ddyAfterLastPeakOrValley = ddy;

% if last peak is before last valley, look for ddy==0 after last valley
if ddyLastPeakInMilliSeconds < ddyLastValleyInMilliSeconds
    ddyAfterLastPeakOrValley(1:ddyLastValleyInDataPoints) = [];
    ddyBasedOffsetInDataPoints = ddyLastValleyInDataPoints + find(round(ddyAfterLastPeakOrValley)==0, 1);
    ddyAfterLastPeakOrValleyInDataPoints = ddyLastValleyInDataPoints;
    
    if isempty(ddyBasedOffsetInDataPoints)
        ddyCrossesZeroPt = find(diff(sign(ddyAfterLastPeakOrValley))>0, 1);
        ddyBasedOffsetInDataPoints = ddyLastValleyInDataPoints + ddyCrossesZeroPt;
    end
        
else
    ddyAfterLastPeakOrValley(1:ddyLastPeakInDataPoints) = [];
    ddyBasedOffsetInDataPoints = ddyLastPeakInDataPoints + find(round(ddyAfterLastPeakOrValley)==0, 1);
    ddyAfterLastPeakOrValleyInDataPoints = ddyLastPeakInDataPoints;
    
    if isempty(ddyBasedOffsetInDataPoints)
        ddyCrossesZeroPt = find(diff(sign(ddyAfterLastPeakOrValley))<0, 1);
        ddyBasedOffsetInDataPoints = ddyLastPeakInDataPoints + ddyCrossesZeroPt;
    end
end 

ddyBasedOffsetInMilliSeconds = xForDdy(ddyBasedOffsetInDataPoints);
ddyAfterLastPeakOrValleyInMilliSeconds = xForDdy(ddyAfterLastPeakOrValleyInDataPoints);


%% Find AP offset based on avgAP==0 after peak
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


%% Calculate AP width
halfWidth = avgAPpeakInMilliSeconds - avgAPvalleyInMilliSeconds;
biphasicDuration = avgAPpeakInMilliSeconds - ddyBasedOnsetInMilliSeconds;
totalDuration = ddyBasedOffsetInMilliSeconds - ddyBasedOnsetInMilliSeconds;
totalDuration2 = avgAPoffsetInMilliSeconds - ddyBasedOnsetInMilliSeconds;


%% Plots
% Plot all APs and avg AP (raw, baseline subtracted)
figure('name', strcat(obj.file,' Baseline Subtracted Non-filtered'));     
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
    title([obj.file ' all APs Raw Baseline Subtracted'],'Interpreter','none');
hold off;

% Plot the first (blue) and second (red) derivative of the avg
figure('name', strcat(obj.file,' y dy ddy Normalized'));
hold on;
    plot(xSubset, avgAP,'Color','black','LineWidth',1);
    plot(xForDy, dy,'Color', 'b', 'LineWidth', 1);
    plot(xForDdy, ddy,'Color', 'r', 'LineWidth', 1);
    plot(ddyAfterLastPeakOrValleyInMilliSeconds, ddy(ddyAfterLastPeakOrValleyInDataPoints), 'o', 'color', 'k');
    plot(ddyBasedOffsetInMilliSeconds, ddy(ddyBasedOffsetInDataPoints), '<', 'color', 'b');
    plot(ddyBasedOnsetInMilliSeconds, ddy(ddyBasedOnsetInDataPoints), '>', 'color', 'b');
    line([0 10],[-ddyValleyThreshold -ddyValleyThreshold], 'LineStyle', '--');
    line([0 10],[ddyPeakThreshold ddyPeakThreshold], 'LineStyle', '--');
    xlabel('Time (ms)');
    ylabel('Amplitude');
    title([obj.file ' Normalized y dy ddy'],'Interpreter','none');
    legend('avg AP', 'avg AP dy', 'avg AP ddy', 'Location', 'northwest');
hold off;
movegui('southeast');

% Plot the first (blue) and second (red) derivative of the avg - Normalized
figure('name', strcat(obj.file,' y dy ddy'));
hold on;
    plot(xSubset, avgAP/max(avgAP),'Color','black','LineWidth',1);
    plot(xForDy, dy/max(dy),'Color', 'b', 'LineWidth', 1);
    plot(xForDdy, ddy/max(ddy),'Color', 'r', 'LineWidth', 1);
    plot(ddyAfterLastPeakOrValleyInMilliSeconds, ddy(ddyAfterLastPeakOrValleyInDataPoints)/max(ddy), 'o', 'color', 'k');
    plot(ddyBasedOffsetInMilliSeconds, ddy(ddyBasedOffsetInDataPoints)/max(ddy), '<', 'color', 'b');
    plot(ddyBasedOnsetInMilliSeconds, ddy(ddyBasedOnsetInDataPoints)/max(ddy), '>', 'color', 'b');
    xlabel('Time (ms)');
    ylabel('Normalized Amplitude (au)');
    title([obj.file ' Normalized y dy ddy'],'Interpreter','none');
    legend('avg AP', 'avg AP dy', 'avg AP ddy', 'Location', 'northwest');
hold off;
movegui('northeast');


%% Store and save data as CSV
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
    totalDuration, ...
    totalDuration2]    

% Save csv file with data 
filename = strcat(obj.file, " - AP width ddy");
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
    'ddyValleyThreshold', ...
    'ddyPeakThreshold', ...
    'nAPtotal', ...
    'ddyBasedOnsetInMilliSeconds', ...
    'avgAPvalleyInMilliSeconds', ...
    'avgAPpeakInMilliSeconds', ...
    'ddyBasedOffsetInMilliSeconds', ...
    'avgAPoffsetInMilliSeconds', ...
    'halfWidth', ...
    'biphasicDuration', ...
    'totalDuration', ...
    'totalDuration2'})
writetable(labeledData,fulldirectory);
disp('I saved the CSV, you save the pics')
   
end    