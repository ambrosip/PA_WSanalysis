% Hello, here's the breakdown of my name/usage:

% lightvsfiring: I am used to analyze and display firing rate before, during and after light stim.

% ON: Use me for cell attached and loose seal recordings.

% bandpass: To clean the data from baseline drift and fast artifacts (such as light stim onset/offset,
% outflow bubbles), I use a bandpass filter.

% AUTO: I am smart code, I identify the light stim parameters
% automatically.

% hist: I will plot APs per sweep for all sweeps and show a histogram of
% firing rate in 1s bins for ALL sweeps. This is to help analyze irregular
% firing cells.

function lightvsfiringONbandpassAUTOhist(obj)
    
%%  USER INPUT

peaksOrValleys = 'v';
highpassThreshold = 100;
lowpassThreshold = 1500;
MinPeakHeight = 9;
ymax = 75;
MinPeakDistance = 0.025;
ymaxhist = 30;
discardedSweeps = 1;

% LightExtensionFactor = 1;
% ZoomWindow = 0.25;
% savefileto = 'D:\CORONAVIRUS DATA\From MATLAB';
    
    % peaksOrValleys: detect peaks (use 'peaks') or valleys ('v'):
    % peaksOrValleys = 'peaks'; 
    % peaksOrValleys = 'v'; 
     
    % highpassThreshold: attenuate low frequency noise below this frequency:
    % highpassThreshold = 100;
    
    % lowpassThreshold: attenuate high frequency noise above this frequency:
    % lowpassThreshold = 1500;
    % lowpassThreshold = 2000;  % used for striatum
    % lowpassThreshold = 1000;  % too low - may cause sinusoidal artifacts
    
    % MinPeakHeight: detect peaks or valleys above/below this amplitude threshold:
    % MinPeakHeight = 8;
    % MinPeakHeight = 5;
    % MinPeakHeight = 20;
    
    % ymax: for illustration purposes, use this value as the max range for the y-axis:
    % ymax = 75;
    
    % MinPeakDistance: to further discard high frequency noise, only look for peaks that are
    % this time interval apart: (note that this interval caps the max
    % signal frequency you can detect)
    % MinPeakDistance = 0.0005;   % used for striatum
    % MinPeakDistance = 0.005;    % used for firing > 50 Hz    
    % MinPeakDistance = 0.05;     % used for firing > 12 Hz 
    
    % ymaxhist: for illustration purposes, use this value as the max range for the
    % y-axis in the histogram (max instantaneous firing frequency):
    % ymaxhist = 15;    
    % ymaxhist = 200;   % used for firing > 100 Hz
    % ymaxhist = 100;   % used for firing > 50 Hz
    % ymaxhist = 50;    % used for firing > 25 Hz
    % ymaxhist = 25;    % used for firing > 12 Hz
    % ymaxhist = 12;    % used for firing < 12 Hz
    
    % LightExtensionFactor: to account for lingering effects after the end of the light pulse,
    % extend the time interval considered under "light effect". To NOT
    % extend the light pulse (my default), set this variable to 1.
    % LightExtensionFactor = 1;
        
    % ZoomWindow is the time before/after LightOnsetTime that will be shown in the zoomed plot
    % ZoomWindow = 0.25;
    
    % savefileto: save CSV files to:
    % savefileto = 'D:\CORONAVIRUS DATA\From MATLAB';
    
    % discardedSweeps: number of sweeps at the end that you want to discard
    % (in case the cell started dying)
    
%% MAIN CODE    
    
    % creating matrixes that will be filled
    allTimeStamps = [];
    
    % getting info from file
    mouseNumber = getMouseNumber(obj);
    experimentDate = getExperimentDate(obj);
    samplingFrequency = obj.header.Acquisition.SampleRate;

    % finding sweep numbers from file name
    [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj); 
    
    % checking for incomplete sweeps and not analyzing incomplete sweeps - to
    % avoid this error: "Index exceeds the number of array elements (0)".      
    if numel(fieldnames(obj.sweeps)) <= obj.header.NSweepsPerRun  
        lastSweepNumber = firstSweepNumber + numel(fieldnames(obj.sweeps)) - discardedSweeps;
        allSweeps = firstSweepNumber:lastSweepNumber;
    end
    
    
%% PLOT FULL HISTOGRAM (0-30 s)
    
    % create figure to illustrate raster plot (APs per time for all sweeps)
    figure('name', strcat(obj.file, ' raster plot'));
    subplot(2,1,1)
    hold on;
    
    % behold this huge 'for loop'
    for sweepNumber = allSweeps
        
        [x,y] = obj.xy(sweepNumber, 1);
        yFiltered = bandpass(y,[highpassThreshold lowpassThreshold],samplingFrequency);
        
        if peaksOrValleys == 'peaks'
            [pks,locs,w,p] = findpeaks(yFiltered,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
        else
            [pks,locs,w,p] = findpeaks(-yFiltered,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
        end
        
        [xch2,ych2] = obj.xy(sweepNumber, 2);  
        
%         sweepDuration = obj.header.Acquisition.Duration;
%         sweepTime=0;
%         inverseISIperSecBin=[];
        
        % finding info about light stim        
        lightPulseStart = find(diff(ych2>1)>0);
        lightPulseEnd = find(diff(ych2<1)>0);
        LightOnsetTime = lightPulseStart(1)/samplingFrequency;                       % in seconds
        StimDur = (lightPulseEnd(end)-lightPulseStart(end))/samplingFrequency;       % duration of each light pulse in the train (s)
        StimInterval = (lightPulseStart(2)-lightPulseStart(1))/samplingFrequency;    % interval between each pulse (s)
        StimFreq = 1/StimInterval;                                                   % frequency of the light stim (Hz)
        LightDur = (lightPulseStart(end)-lightPulseStart(1))/samplingFrequency + StimInterval;  % duration of the whole light train stim (s)
        
        % create list of sweepNumber with the same size as list of timestamps for plotting
        sweepNumberArray = sweepNumber.* ones(length(locs),1);
        
        % plotting raster plot 
        plot(locs, sweepNumberArray, '|', 'Color', 'k');
        axis([0 30 firstSweepNumber-1 lastSweepNumber+1])
        ylabel(strcat('Sweeps (', num2str(length(allSweeps)), ')'));
        yticks([]);
        xticks([]);

        % getting all of the APs timestamps
        allTimeStamps = [allTimeStamps; locs]

    end
    
    % adding light stim
    rectangle('Position', [LightOnsetTime firstSweepNumber-1 LightDur lastSweepNumber+1], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');
    title([strcat(obj.file, ' raster plot')],'Interpreter','none');  
    
    % stop plotting things on this subplot
    hold off;
    
    % flip the y-axis so that the first sweep is at the top and the last
    % sweep is at the bottom
    set(gca, 'YDir','reverse');

    % counting APs accross all sweeps
    edges = [0:30];
    [N, edges] = histcounts(allTimeStamps,edges);
    firingHz = N/length(allSweeps);
        
    % plotting Hz histogram 
    subplot(2,1,2)
    hold on;
    histogram('BinEdges',0:30,'BinCounts',firingHz,'EdgeColor','none','FaceColor','black','FaceAlpha',1);
    rectangle('Position', [LightOnsetTime 0 LightDur ymaxhist], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');
    xlabel('Time (s)');
    ylabel('Firing rate (Hz)');
    axis([0 30 0 ymaxhist])
    set(gcf,'Position',[200 200 500 400]);
    xticks([0 30]);
    yticks([0 ymaxhist]);
    hold off;
    
    
%% PLOT ZOOMED IN HISTOGRAM (14-21 s)
    
    % create figure to illustrate raster plot (APs per time for all sweeps)
    figure('name', strcat(obj.file, ' raster plot zoom'));
    subplot(2,1,1)
    hold on;
    
    % going through all sweeps
    for sweepNumber = allSweeps
        
        [x,y] = obj.xy(sweepNumber, 1);
        yFiltered = bandpass(y,[highpassThreshold lowpassThreshold],samplingFrequency);
        
        if peaksOrValleys == 'peaks'
            [pks,locs,w,p] = findpeaks(yFiltered,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
        else
            [pks,locs,w,p] = findpeaks(-yFiltered,x,'MinPeakHeight',MinPeakHeight,'MinPeakDistance',MinPeakDistance);
        end
        
        [xch2,ych2] = obj.xy(sweepNumber, 2);  
               
        % create list of sweepNumber with the same size as list of timestamps for plotting
        sweepNumberArray = sweepNumber.* ones(length(locs),1);
        
        % plotting raster plot 
        plot(locs, sweepNumberArray, '|', 'Color', 'k');
        axis([LightOnsetTime-LightDur LightOnsetTime+2*LightDur firstSweepNumber-1 lastSweepNumber+1])
        ylabel(strcat('Sweeps (', num2str(length(allSweeps)), ')'));
        yticks([]);
        xticks([]);

    end
    
    % adding light stim
    rectangle('Position', [LightOnsetTime firstSweepNumber-1 LightDur lastSweepNumber+1], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');
    title([strcat(obj.file, ' raster plot zoom')],'Interpreter','none');  
    
    % stop plotting things on this subplot
    hold off;
    
    % flip the y-axis so that the first sweep is at the top and the last
    % sweep is at the bottom
    set(gca, 'YDir','reverse');

    % counting APs accross all sweeps
    edges = [0:30];
    [N, edges] = histcounts(allTimeStamps,edges);
    firingHz = N/length(allSweeps);
        
    % plotting Hz histogram 
    subplot(2,1,2)
    hold on;
    histogram('BinEdges',0:30,'BinCounts',firingHz,'EdgeColor','none','FaceColor','black','FaceAlpha',1);
    rectangle('Position', [LightOnsetTime 0 LightDur ymaxhist], 'FaceColor', [0 0.4470 0.7410 0.1], 'EdgeColor', 'none');
    xlabel('Time (s)');
    ylabel('Firing rate (Hz)');
    axis([LightOnsetTime-LightDur LightOnsetTime+2*LightDur 0 ymaxhist])
    set(gcf,'Position',[200 200 500 400]);
%     xticks([0 30]);
    yticks([0 ymaxhist]);
    hold off;
    
end