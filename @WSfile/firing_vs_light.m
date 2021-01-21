%{ 
DOCUMENTATION
Created: 2021 01 19
Author: PA

This function is used to analyze light-evoked changes in firing rate.
This function automatically identifies the light stim parameters.

INPUTS explained:
    - discardedSweeps: number of sweeps from the end that will not be
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

    - ZoomWindow: time before/after LightOnsetTime that will be shown in
    the zoomed plot.

    - ymaxIsiCV: for illustration purposes, use this value as the max range
    for the y-axis in the ISI histogram.

    - savefileto: save CSV files to this folder.

INPUTS defaults:
    discardedSweeps = 1;
    peaksOrValleys = 'v';   
    highpassThreshold = 100;
    lowpassThreshold = 1500;    
    minPeakHeight = 15;         
    minPeakDistance = 0.025;    
    LightExtensionFactor = 1;
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

TO DO:
    - add all plots - sweep by sweep plots, AP width
    - complete documentation

%}

function firing_vs_light(obj)

%%  USER INPUT

% Affects data analysis:
discardedSweeps = 1;
peaksOrValleys = 'v';   
highpassThreshold = 100;
lowpassThreshold = 1500;    
minPeakHeight = 15;         
minPeakDistance = 0.025;    
lightExtensionFactor = 1;

% Affects data display:
ymax = 75;
ymaxhist = 15;
zoomWindow = 0.25;
ymaxIsiCV = 150;

% Affects data saving:
savefileto = 'D:\CORONAVIRUS DATA\From MATLAB';


%% PREP - get info from file and create arrays

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% getting info from file
mouseNumber = getMouseNumber(obj);
experimentDate = getExperimentDate(obj);
samplingFrequency = obj.header.Acquisition.SampleRate;

% getting sweep numbers from file name
[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj); 

% checking for incomplete sweeps and not analyzing incomplete sweeps - to
% avoid this error: "Index exceeds the number of array elements (0)".      
if numel(fieldnames(obj.sweeps)) <= obj.header.NSweepsPerRun  
    lastSweepNumber = firstSweepNumber + numel(fieldnames(obj.sweeps)) - discardedSweeps;
    allSweeps = firstSweepNumber:lastSweepNumber;
end

% creating matrixes/arrays that will be filled
allTimeStamps = [];
baselineTimeStamps = [];
allIsiBaseline = [];
hzBySweep = [];
isiCvBySweep = [];
data = [];

% creating cell arrays that will be filled (can take columns of different sizes)
tsBySweep = {};
isiBySweep = {};
sweepNumberArrayBySweep = {};


%% SWEEP BY SWEEP ANALYSIS

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


    % getting peak/valley amplitudes pre, during and post light
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

    
    % getting timestamps pre, during, and post light
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
    

    % storing sweep by sweep data in a cell array
    % to export later
    % fyi access sweep 1 data using cell2mat(sweepData(1))
    tsBySweep = [tsBySweep, locs];
    isiBySweep = [isiBySweep, isiBaseline];
    sweepNumberArrayBySweep = [sweepNumberArrayBySweep, sweepNumberArray];
    
    % storing sweep by sweep data in an array for easy mean & std calculations later
    hzBySweep = [hzBySweep, length(locsBaseline)/lightOnsetTime];
    isiCvBySweep = [isiCvBySweep, std(isiBaseline)/mean(isiBaseline)];
    
    % storing a subset of sweep by sweep data that will be exported
    % this dataset is the one I'm used to manipulating with a few changes
    data = [data; ...
        mouseNumber, ...
        experimentDate, ...
        sweepNumber, ...
        discardedSweeps, ...
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
        length(locsPreLight), ...
        length(locsDuringLight), ...
        length(locsPostLight), ...
        length(locsPreLight)/lightDur, ...
        length(locsDuringLight)/lightDur, ...
        length(locsPostLight)/lightDur];    
    
end


%% CELL ANALYSIS

hzMean = mean(hzBySweep);
hzStd = std(hzBySweep);

% counting APs accross all sweeps
edges = [0:30];
[N, edges] = histcounts(allTimeStamps,edges);
firingHz = N/length(allSweeps);


%% PLOT 1 - ISI CV

% ISI CV pre-light across all sweeps
figure('name', strcat(obj.file, " ", analysisDate, ' - firing_vs_light - baseline ISI CV'));

% edges: from 0-1s in 1ms steps
edges = [0:0.001:1]; 

histogram(allIsiBaseline, edges);
title([strcat(obj.file, ' baseline ISI CV')],'Interpreter','none');
xlabel('ISI (s)');
ylabel('Counts per bin');
axis([0 1 0 ymaxIsiCV])
xticks([0 1]);
set(gcf,'Position',[50 550 400 400]);
yticks([0 ymaxIsiCV]);


%% PLOT 1.2 - ISI CV free y range

% ISI CV pre-light across all sweeps
figure('name', strcat(obj.file, " ", analysisDate, ' - firing_vs_light - baseline ISI CV yrange'));

% edges: from 0-1s in 1ms steps
edges = [0:0.001:1]; 

histogram(allIsiBaseline, edges);
title([strcat(obj.file, ' baseline ISI CV yrange')],'Interpreter','none');
xlabel('ISI (s)');
ylabel('Counts per bin');
axis([0 0.5 0 inf])
xticks([0 0.5]);
set(gcf,'Position',[50 50 400 400]);


%% PLOT 2 - Raster plot and histogram - Zoomed in

figure('name', strcat(obj.file, " ", analysisDate, ' - firing_vs_light - raster and hist zoom'));
subplot(2,1,1)
hold on;

% plotting raster plot 
for sweep = 1:length(allSweeps)
    plot(cell2mat(tsBySweep(sweep)), cell2mat(sweepNumberArrayBySweep(sweep)), '|', 'Color', 'k')
end

% zooming in and beautifying raster plot
title([strcat(obj.file, ' raster and hist zoom')],'Interpreter','none');
axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur firstSweepNumber-1 lastSweepNumber+1])
ylabel(strcat('Sweeps (', num2str(length(allSweeps)), ')'));
yticks([]);
xticks([]);

% adding light stim
rectangle('Position', [lightOnsetTime firstSweepNumber-1 lightDur lastSweepNumber+1], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');

% stop plotting things on this subplot
hold off;

% flip the y-axis so that the first sweep is at the top and the last
% sweep is at the bottom
set(gca, 'YDir','reverse');

% plotting Hz histogram 
subplot(2,1,2)
hold on;

histogram('BinEdges', 0:30, 'BinCounts', firingHz, 'DisplayStyle', 'stairs', 'EdgeColor', 'k'); % 'EdgeColor','none',

% plot light stim as rectangle
rectangle('Position', [lightOnsetTime 0 lightDur ymaxhist], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');

% plot Hz mean as horizontal line
yline(hzMean, '--');

% plot +- 2 SD as rectangle around mean
% [x y width height]
rectangle('Position', [0 hzMean-(2*hzStd) 30 4*hzStd], 'FaceColor', [0 0 0 0.1], 'EdgeColor', 'none');

xlabel('Time (s)');
ylabel('Firing rate (Hz)');
axis([lightOnsetTime-lightDur lightOnsetTime+2*lightDur 0 ymaxhist])
set(gcf,'Position',[500 50 500 400]);
%     xticks([0 30]);
yticks([0 ymaxhist]);
hold off;
% movegui('east')


%% PLOT 3 - Raster plot and histogram - Complete

figure('name', strcat(obj.file, " ", analysisDate, ' - firing_vs_light - raster and hist'));
subplot(2,1,1)
hold on;

% plotting raster plot 
for sweep = 1:length(allSweeps)
    plot(cell2mat(tsBySweep(sweep)), cell2mat(sweepNumberArrayBySweep(sweep)), '|', 'Color', 'k')
end

% beautifying raster plot
title([strcat(obj.file, ' raster and hist')],'Interpreter','none');
axis([0 30 firstSweepNumber-1 lastSweepNumber+1])
ylabel(strcat('Sweeps (', num2str(length(allSweeps)), ')'));
yticks([]);
xticks([]);

% adding light stim
rectangle('Position', [lightOnsetTime firstSweepNumber-1 lightDur lastSweepNumber+1], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');

% stop plotting things on this subplot
hold off;

% flip the y-axis so that the first sweep is at the top and the last
% sweep is at the bottom
set(gca, 'YDir','reverse');

% plotting Hz histogram 
subplot(2,1,2)
hold on;

histogram('BinEdges', 0:30, 'BinCounts', firingHz, 'DisplayStyle', 'stairs', 'EdgeColor', 'k'); % 'EdgeColor','none',

% plot light stim as rectangle
rectangle('Position', [lightOnsetTime 0 lightDur ymaxhist], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');

% plot Hz mean as horizontal line
yline(hzMean, '--');

% plot +- 2 SD as rectangle around mean
% [x y width height]
rectangle('Position', [0 hzMean-(2*hzStd) 30 4*hzStd], 'FaceColor', [0 0 0 0.1], 'EdgeColor', 'none');

xlabel('Time (s)');
ylabel('Firing rate (Hz)');
axis([0 30 0 ymaxhist])
set(gcf,'Position',[500 550 500 400]);
% xticks([0 30]);
yticks([0 ymaxhist]);
hold off;



%% EXPORTING XLS files

% stores key sweep by sweep data
filename = strcat(obj.file, " ", analysisDate, " - firing_vs_light - sweep_by_sweep");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(data);
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...
    {'mouse', ...
    'date', ...
    'sweep', ...
    'discardedSweeps', ...
    'peaks(1)OrValleys(-1)', ...
    'highpassThreshold', ...
    'lowpassThreshold', ...
    'minPeakHeight', ...        
    'minPeakDistance', ...    
    'lightExtensionFactor', ...
    'lightPulseDur(s)', ... 
    'lightStimFreq(Hz)', ...
    'lightDur(s)', ...
    'baselineHz', ...
    'baselineIsiMean(s)', ...
    'baselineIsiSD(s)', ...
    'baselineIsiCV', ...
    'preLightAPs', ...
    'duringLightAPs', ...
    'postLightAPs', ...
    'preLightHz', ...
    'duringLightHz', ...
    'postLightHz'});
writetable(labeledData,fulldirectory);
disp('I saved the sweep_by_sweep xls file')

% stores all timestamps from all sweeps in a single column
filename = strcat(obj.file, " ", analysisDate, " - firing_vs_light - all_AP_timestamps");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
writematrix(allTimeStamps,fulldirectory);
disp('I saved the all_AP_timestamps xls file')


end








    
    
    
    
    
    
    
    
    
    


