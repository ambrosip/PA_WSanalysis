%{ 
DOCUMENTATION
Created: 2021 09 19
Author: PA

This function is used to analyze the kinetics of post-synaptic currents
(PSCs), either IPSCs, EPSCs, or ChR2-mediated currents.

This code was adapted from psc_vs_light

INPUTS explained:
    TO DO

INPUTS defaults:
    % Affects data analysis:
    lightStimCh = 2;
    discardedSweeps = [];
    discardedSweepsFromEnd = 1;
    inwardORoutward = 1;    % 1 (positive) is outward; -1 (negative) in inward
    baselineDurationInSeconds = 0.5;
    lightPulseAnalysisWindowInSeconds = 0.02;
    thresholdInDataPts = 10;
    rsTestPulseOnsetTime = 1;

    % Affects data display:
    ymax = 150;

    % Affects data saving:
    savefileto = 'D:\CORONAVIRUS DATA\From MATLAB';
    
OUTPUTS:
    TO DO

ASSUMPTIONS: 
    TO DO

TO DO:
    TO DO
%}

function psc_vs_light_single(obj)
%%%  USER INPUT ==================================================

% Affects data analysis:
lightStimCh = 2;
discardedSweeps = [];
discardedSweepsFromEnd = 0;
inwardORoutward = -1;    % 1 (positive) is outward; -1 (negative) in inward
baselineDurationInSeconds = 0.01;
lightPulseAnalysisWindowInSeconds = 0.02;   %% ALERT! Changed from 0.02 to 0.015 on 2023/3/4   
thresholdInDataPts = 8; %% ALERT! Changed from 10 to 5; now changed to 8 on 2023/3/4 
rsTestPulseOnsetTime = 1; %% ALERT! Changed from 1 to 0.1

% Affects decay fit
bGuess = 1;         % -1000     % 1
fastTauGuess = 0.01;    % in seconds
slowTauGuess = 0.1;     % in seconds

% Affects data display:
ymin = -12000;      % -4200   % -3600   % -2050   % -375     % -1500
ymax = 600;        % 700     % 600     % 50      % 150      % 600

% Affects data saving:
savefileto = 'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries\2023\2023 08 24 scracm des m117';


%% PREP - get info from file and create arrays ==================

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% getting info from file
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

% calculating variables based on user input
baselineDurationInDataPts = baselineDurationInSeconds * samplingFrequency;
lightPulseAnalysisWindowInDataPts = lightPulseAnalysisWindowInSeconds * samplingFrequency;
rsBaselineDataPointInterval = ((rsTestPulseOnsetTime-0.05)*samplingFrequency):(rsTestPulseOnsetTime*samplingFrequency);
rsFirstTransientDataPointInterval = (rsTestPulseOnsetTime*samplingFrequency):(rsTestPulseOnsetTime+0.0025)*samplingFrequency;
mouseNumber = getMouseNumber(obj);
experimentDate = getExperimentDate(obj);

% creating matrixes
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


%% SWEEP BY SWEEP ANALYSIS ======================================

% get data from all sweeps in file
for sweepNumber = allSweeps
    
    % get light stim data
    [xch2,ych2] = obj.xy(sweepNumber, lightStimCh);      
    
    % get light stim parameters
    % I commented out the stuff cuz this code applies to a single light pulse
    lightPulseStart = find(diff(ych2>1)>0);                                      % list of light pulse onset data points
    lightPulseEnd = find(diff(ych2<1)>0);                                        % list of light pulse offset data points
    lightOnsetTime = lightPulseStart(1)/samplingFrequency;                       % onset of first light pulse in seconds
    stimDur = (lightPulseEnd(end)-lightPulseStart(end))/samplingFrequency       % duration of each light pulse in the train (s)
    % commented out stuff below because it does not apply to a single light
    % stim:
%     stimInterval = (lightPulseStart(2)-lightPulseStart(1))/samplingFrequency;    % interval between each pulse (s)
%     stimFreq = 1/stimInterval;                                                   % frequency of the light stim (Hz)
%     lightDur = (lightPulseStart(end)-lightPulseStart(1))/samplingFrequency + stimInterval;  % duration of the whole light train stim (s)  
    
    xmin = lightOnsetTime-0.02;
    xmax = lightOnsetTime+0.2;

    % get raw data
    [x,y] = obj.xy(sweepNumber, 1);
    
    % baseline subtraction
    baselineStart = lightPulseStart(1) - baselineDurationInDataPts;
    yBaselineSub = y-mean(y(baselineStart:lightPulseStart(1)));
    baselineCurrent = mean(y(baselineStart:lightPulseStart(1)));
    
    % saving data for niceplot
    % y data for each sweep is in a column
    yBaselineSubAll = [yBaselineSubAll, yBaselineSub]; 
    baselineCurrentAll = [baselineCurrentAll, baselineCurrent];
    
%     % filter data
%     yFiltered = bandpass(y,[highpassThreshold lowpassThreshold],samplingFrequency);
%     
%     % find peaks or valleys based on user input
%     if peaksOrValleys == 'peaks'
%         [pks,locs,w,p] = findpeaks(yFiltered,x,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);
%     else
%         [pks,locs,w,p] = findpeaks(-yFiltered,x,'MinPeakHeight',minPeakHeight,'MinPeakDistance',minPeakDistance);
%     end

    % calculating series resistance
    rsBaselineCurrent = mean(y(rsBaselineDataPointInterval));
    rsTransientCurrent = min(y(rsFirstTransientDataPointInterval));
    dCurrent = rsTransientCurrent-rsBaselineCurrent;
    dVoltage = -5;
    seriesResistance = 1000*dVoltage/dCurrent; %mV/pA equals Gohm
    
    % put series resistance values from each sweep into a different column
    allRs = [allRs, seriesResistance];
    
    % clean up matrices for that will be used in the next loop
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
    
    % Data that will be exported (each sweep is a row)       
    data = [data; ...
        mouseNumber, ...
        experimentDate, ...
        sweepNumber, ...
        discardedSweepsFromEnd, ...
        inwardORoutward, ...
        baselineDurationInSeconds, ...
        lightPulseAnalysisWindowInSeconds, ...
        thresholdInDataPts, ...
        rsTestPulseOnsetTime, ...
        stimDur, ... 
        seriesResistance, ...
        baselineCurrent];     
        
end

% add extra columns for storage
data = [data, ... 
    lightEvokedCurrentsAllSweeps, ...
    allLightEvokedResponseOnsetLatencyInMilliSecondsAllSweeps, ...
    allLightEvokedResponsePeakLatencyInMilliSecondsAllSweeps, ...
    allTimeTo10percentOfPeakInMilliSecondsAllSweeps, ...
    allTimeTo90percentOfPeakInMilliSecondsAllSweeps, ...
    allRiseTimeInMilliSecondsAllSweeps];

% WARNING: I turned this light detection-based x axis into input-defined
% xmin = lightOnsetTime - baselineDurationInSeconds;
% xmax = lightOnsetTime + baselineDurationInSeconds;
% % To accommodate single light pulse, I removed lightDur from the xmax calculation!
% % xmax = lightOnsetTime + lightDur + baselineDurationInSeconds;


%% CELL Analysis

yBaselineSubAllMean = mean(yBaselineSubAll,2);

% ASSUMPTION ALERT
% Assuming that all sweeps have the same, single o-stim
pulseOnset = lightPulseStart(1)
% pulseOnset = lightPulseStart.'
afterLightDataPoint = pulseOnset + lightPulseAnalysisWindowInDataPts;

% get amplitude and location (index) of lightEvokedCurrents
% for rise time calculation, find index of 10% and 90% of peak
% look for a threshold cross instead of an exact match - find(a<=b) instead of find(a==b) 
% outward current is +1
% inward current is 0 or -1
if inwardORoutward == 1 
    [lightEvokedCurrentAmp, lightEvokedCurrentLoc] = max(yBaselineSubAllMean(pulseOnset:afterLightDataPoint))
    latencyTo10percentOfPeakInDataPoints = find(yBaselineSubAllMean(pulseOnset:afterLightDataPoint) >= 0.1 * lightEvokedCurrentAmp, 1) + pulseOnset;
    latencyTo90percentOfPeakInDataPoints = find(yBaselineSubAllMean(pulseOnset:afterLightDataPoint) >= 0.9 * lightEvokedCurrentAmp, 1) + pulseOnset;
else
    [lightEvokedCurrentAmp, lightEvokedCurrentLoc] = min(yBaselineSubAllMean(pulseOnset:afterLightDataPoint))
    latencyTo10percentOfPeakInDataPoints = find(yBaselineSubAllMean(pulseOnset:afterLightDataPoint) <= 0.1 * lightEvokedCurrentAmp, 1) + pulseOnset;
    latencyTo90percentOfPeakInDataPoints = find(yBaselineSubAllMean(pulseOnset:afterLightDataPoint) <= 0.9 * lightEvokedCurrentAmp, 1) + pulseOnset;
end

% adjust latency to peak Loc
lightEvokedCurrentLoc = lightEvokedCurrentLoc + pulseOnset

% find onset latency
onsetLatencyInDataPoints = min(find(conv2(sign(diff(yBaselineSubAllMean(pulseOnset:afterLightDataPoint))), zeros(thresholdInDataPts,1)+(inwardORoutward), 'valid')==thresholdInDataPts));
onsetLatencyInDataPoints = onsetLatencyInDataPoints + pulseOnset;

% for decay tau calculation, find index of return to baseline after peak
latencyToZeroAfterPeakInDataPoints = find(round(yBaselineSubAllMean(lightEvokedCurrentLoc:end)) >= 0, 1) + lightEvokedCurrentLoc

% preparing to fit double-term exponential model to the decay data - create y subset
xminFit = lightEvokedCurrentLoc;
xmaxFit = latencyToZeroAfterPeakInDataPoints;
yBaselineSubAllMeanSubset = yBaselineSubAllMean;
yBaselineSubAllMeanSubset(latencyToZeroAfterPeakInDataPoints:end) = [];
yBaselineSubAllMeanSubset(1:lightEvokedCurrentLoc) = [];

% preparing to fit double-term exponential model to the decay data - create x subset
xSubset = x;
xSubset(latencyToZeroAfterPeakInDataPoints:end) = [];
xSubset(1:lightEvokedCurrentLoc) = [];
xSubset = round(xSubset,6) - round(lightEvokedCurrentLoc/samplingFrequency,6);

% selecting a double exponential for decay fitting
g = fittype('a*exp(-x/fastTau) + b*exp(-x/slowTau)');

% fit double-term exponential model to the decay data 
% If you don't give MATLAB a StartPoint, it fails horribly
% StartPoint order: [a, b, fastTau, slowTau]
% You can check the order by using this line:
% coefficientNames = coeffnames(g)
% I guessed StartPoint values by manually adjusting a double exponential
% The single exponential was NOT giving me a good fit
[decayFit,gof,output] = fit(xSubset,yBaselineSubAllMeanSubset,g,'StartPoint',[lightEvokedCurrentAmp, bGuess, fastTauGuess, slowTauGuess]);

% save coefficients a and b
decayFitA = decayFit.a;
decayFitB = decayFit.b;

% % For troubleshooting and estimating StartPoints:
% figure;
% plot(decayFit);
% figure;
% hold on;
% plot(xSubset+round(lightEvokedCurrentLoc/samplingFrequency,6), decayFit(xSubset))
% plot(x, yBaselineSubAllMean);
% plot(xSubset+round(lightEvokedCurrentLoc/samplingFrequency,6), lightEvokedCurrentAmp*exp(-xSubset/0.05));
% hold off;

% convert data points to ms
onsetLatencyInMilliSeconds = 1000 * (onsetLatencyInDataPoints - pulseOnset) / samplingFrequency
latencyToPeakInMilliSeconds = 1000 * (lightEvokedCurrentLoc - pulseOnset) / samplingFrequency;
latencyToZeroAfterPeakInMilliSeconds = 1000 * (latencyToZeroAfterPeakInDataPoints - pulseOnset) / samplingFrequency;
latencyTo10percentOfPeakInMilliSeconds = 1000 * (latencyTo10percentOfPeakInDataPoints - pulseOnset) / samplingFrequency;
latencyTo90percentOfPeakInMilliSeconds = 1000 * (latencyTo90percentOfPeakInDataPoints - pulseOnset) / samplingFrequency;

% check if onset latency could be calculated to avoid errors
if isempty(onsetLatencyInDataPoints)
    onsetLatencyInMilliSeconds = NaN;
    onsetLatencyInDataPoints = NaN;
end

% convert s to ms
decayFitFastTau = decayFit.fastTau * 1000
decayFitSlowTau = decayFit.slowTau * 1000

% calculate riseTime in MilliSeconds
riseTimeInMilliSeconds = latencyTo90percentOfPeakInMilliSeconds - latencyTo10percentOfPeakInMilliSeconds

% store all data
dataCell = [mouseNumber, ...
    experimentDate, ...
    firstSweepNumber, ...
    lastSweepNumber, ...
    length(allSweeps), ...
    bGuess, ...
    fastTauGuess, ...
    slowTauGuess, ...
    mean(allRs), ...
    min(allRs), ...
    max(allRs), ...
    mean(baselineCurrentAll), ...
    min(baselineCurrentAll), ...
    max(baselineCurrentAll), ...
    mean(lightEvokedCurrentsAllSweeps), ...
    std(lightEvokedCurrentsAllSweeps), ...
    mean(allLightEvokedResponseOnsetLatencyInMilliSecondsAllSweeps, 'omitnan'), ...
    std(allLightEvokedResponseOnsetLatencyInMilliSecondsAllSweeps, 'omitnan'), ...
    mean(allLightEvokedResponsePeakLatencyInMilliSecondsAllSweeps, 'omitnan'), ...
    std(allLightEvokedResponsePeakLatencyInMilliSecondsAllSweeps, 'omitnan'), ...    
    mean(allTimeTo10percentOfPeakInMilliSecondsAllSweeps, 'omitnan'), ...
    std(allTimeTo10percentOfPeakInMilliSecondsAllSweeps, 'omitnan'), ...
    mean(allTimeTo90percentOfPeakInMilliSecondsAllSweeps, 'omitnan'), ...
    std(allTimeTo90percentOfPeakInMilliSecondsAllSweeps, 'omitnan'), ...
    mean(allRiseTimeInMilliSecondsAllSweeps, 'omitnan'), ...
    std(allRiseTimeInMilliSecondsAllSweeps, 'omitnan'), ...
    lightEvokedCurrentAmp, ...
    onsetLatencyInMilliSeconds, ... 
    latencyToPeakInMilliSeconds, ... 
    latencyTo10percentOfPeakInMilliSeconds, ... 
    latencyTo90percentOfPeakInMilliSeconds, ... 
    latencyToZeroAfterPeakInMilliSeconds, ... 
    riseTimeInMilliSeconds, ...  
    decayFitFastTau, ...
    decayFitSlowTau, ...
    decayFitA, ...
    decayFitB, ...
    gof.rsquare];


%% Color-blind diverging color scheme

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
                       

%% PLOT - light-evoked currents amplitude ==========================
% 
% f1=figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_single - oPSCs amplitude all')); % naming figure file
% hold on;
% colororder(f1,brownToPurpleHex);
% plot(lightEvokedCurrentsAllSweeps.');
% colormap(brownToPurpleRgb(1:length(allSweeps),:))
% c = colorbar('Ticks',[1,length(allSweeps)]);
% caxis([1 length(allSweeps)]);
% c.Label.String = 'Sweep #';
% set(c, 'YDir', 'reverse');
% line([0 60],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
% axis([-inf inf -ymax ymax]);
% title([fileName ' - psc_vs_light_single - oPSC amplitude all'],'Interpreter','none');
% ylabel(strcat(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
% xlabel('Light pulse');
% hold off;
% movegui('northwest');


%% PLOT - light-evoked currents latency ===========================
% 
% f2=figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_single - oPSCs latency all')); % naming figure file
% hold on;
% colororder(f2,brownToPurpleHex);
% plot(allLightEvokedResponseOnsetLatencyInMilliSecondsAllSweeps.','-o');
% colormap(brownToPurpleRgb(1:length(allSweeps),:))
% c = colorbar('Ticks',[1,length(allSweeps)]);
% caxis([1 length(allSweeps)]);
% c.Label.String = 'Sweep #';
% set(c, 'YDir', 'reverse');
% line([0 60],[5, 5],'Color','black','LineStyle','--')
% line([0 60],[1, 1],'Color','red','LineStyle','--')
% axis([-inf inf 0 10]);
% title([fileName ' - psc_vs_light_single - oPSC latency all'],'Interpreter','none');
% ylabel('Response Latency (ms)');
% xlabel('Light pulse');
% hold off;
% movegui('north');


%% PLOT - rs =====================================================

figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_single - Rs all')); % naming figure file
plot(allSweeps, allRs,'-o');
% plot lines marking 30% increase and 30% decrese in Rs compared to first
% test pulse
line([allSweeps(1) allSweeps(end)],[allRs(1)*0.7, allRs(1)*0.7],'Color','black','LineStyle','--')
line([allSweeps(1) allSweeps(end)],[allRs(1)*1.3, allRs(1)*1.3],'Color','black','LineStyle','--')
axis([allSweeps(1) inf 0 60])
ylabel('Rs (M\Omega)');
xlabel('Sweeps');
title([obj.file ' rs'],'Interpreter','none');
movegui('northeast');


%% PLOT - niceplot of all sweeps - COLORS ================================

% plotting niceplot     
figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_single - Niceplot all colors'));
hold on;
for sweep = 1:size(yBaselineSubAll,2)
    plot(x, yBaselineSubAll(:,sweep),'Color',brownToPurpleRgbTransparent(sweep,:));
end
colormap(brownToPurpleRgb(1:sweep,:))
c = colorbar('Ticks',[1,sweep]);
caxis([1 sweep]);
c.Label.String = 'Sweep #';
set(c, 'YDir', 'reverse');
% plot(x, mean(yBaselineSubAll,2),'Color','black','LineWidth',1.5); 
line([xmin xmax],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
axis([xmin xmax ymin ymax]);
xlabel('Time (s)');
ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
title([fileName ' - psc_vs_light_single - niceplot all colors'],'Interpreter','none');
set(gca,'Visible','off');
set(gcf,'Position',[1400 550 500 400]);

% adding light stim - individual pulses   
% note that this code will use light stim parameters from the last sweep!
% if light stim is not the same accross all sweeps, this will be
% misleading!
for nStim=1:length(lightPulseStart)
    line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-inwardORoutward*(ymax),-inwardORoutward*(ymax)],'Color',[0 0.4470 0.7410],'LineWidth',10)
end

% adding scale bar
xmaxScale = xmax;
xminScale = xmin;
line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
text(xmaxScale-(xmaxScale-xminScale)/8,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
text(xmaxScale-(xmaxScale-xminScale)/6,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/7)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
hold off;
movegui('southeast');


%% PLOT - niceplot of all sweeps and/or MEAN - BLACK ================================

% plotting niceplot     
figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_single - niceplot all'));
hold on;
plot(x, yBaselineSubAll,'Color',[0, 0, 0, 0.25]);
plot(x, mean(yBaselineSubAll,2),'Color','black','LineWidth',0.7); 
line([xmin xmax],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
axis([xmin xmax ymin ymax]);
xlabel('Time (s)');
ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
title([fileName ' - psc_vs_light_single - niceplot all'],'Interpreter','none');
set(gca,'Visible','off');
set(gcf,'Position',[1400 550 500 400]);

% adding light stim - individual pulses   
% note that this code will use light stim parameters from the last sweep!
% if light stim is not the same accross all sweeps, this will be
% misleading!
for nStim=1:length(lightPulseStart)
    line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-inwardORoutward*(ymax),-inwardORoutward*(ymax)],'Color',[0 0.4470 0.7410],'LineWidth',10)
end

% adding scale bar
xmaxScale = xmax;
xminScale = xmin;
line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
text(xmaxScale-(xmaxScale-xminScale)/9,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
text(xmaxScale-(xmaxScale-xminScale)/7,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/7)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
hold off;
movegui('south');


%% PLOT - subtracted baseline current =====================================================

figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_single - baseline current')); % naming figure file
plot(allSweeps, baselineCurrentAll,'-o');
axis([allSweeps(1) inf -300 300])
ylabel('Subtracted Baseline Current (pA)');
xlabel('Sweeps');
title([obj.file ' rs'],'Interpreter','none');
movegui('southwest');


%% PLOT - oPSC amplitude mean +- SD ================================

% figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_single - mean and SD'));
% hold on;
% errorbar(mean(lightEvokedCurrentsAllSweeps), std(lightEvokedCurrentsAllSweeps), 'color', [0 0 0], 'CapSize', 0);
% % errorbar(mean(lightEvokedCurrentsAllSweeps), std(lightEvokedCurrentsAllSweeps),'-o', 'MarkerFaceColor', 'k', 'color', [0 0 0], 'CapSize', 0);
% line([0 60],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
% axis([-inf inf ymin ymax]);
% ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
% title([fileName ' - psc_vs_light_single - mean and SD'],'Interpreter','none');
% xlabel('Light pulse');
% hold off;
% movegui('northwest');


%% PLOT - niceplot with kinetics - BLACK ================================

% plotting niceplot     
figure('name', strcat(fileName, " ", analysisDate, ' - psc_vs_light_single - Niceplot kinetics black'));
hold on;
% plot(x, yBaselineSubAll,'Color',[0, 0, 0, 0.25]);
plot(x, mean(yBaselineSubAll,2),'Color','black','LineWidth',0.7); 
plot(xSubset+round(lightEvokedCurrentLoc/samplingFrequency,6), decayFit(xSubset))
plot(lightEvokedCurrentLoc/samplingFrequency, lightEvokedCurrentAmp, 'o', 'Color', 'red');
line([onsetLatencyInDataPoints/samplingFrequency onsetLatencyInDataPoints/samplingFrequency], [ymin ymax], 'Color',[0.5 0.5 0.5],'LineStyle','--');
line([latencyTo10percentOfPeakInDataPoints/samplingFrequency latencyTo10percentOfPeakInDataPoints/samplingFrequency], [ymin ymax], 'Color',[0.5 0.5 0.5],'LineStyle','--');
line([latencyTo90percentOfPeakInDataPoints/samplingFrequency latencyTo90percentOfPeakInDataPoints/samplingFrequency], [ymin ymax], 'Color',[0.5 0.5 0.5],'LineStyle','--');
line([latencyToZeroAfterPeakInDataPoints latencyToZeroAfterPeakInDataPoints], [ymin ymax], 'Color',[0.5 0.5 0.5],'LineStyle','--');
line([xmin xmax],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
axis([xmin xmax ymin ymax]);
xlabel('Time (s)');
ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
title([fileName ' - psc_vs_light_single - niceplot kinetics'],'Interpreter','none');
set(gca,'Visible','off');
set(gcf,'Position',[1400 550 500 400]);

% adding light stim - individual pulses   
% note that this code will use light stim parameters from the last sweep!
% if light stim is not the same accross all sweeps, this will be
% misleading!
for nStim=1:length(lightPulseStart)
    line([(lightPulseStart(nStim)/samplingFrequency),(lightPulseStart(nStim)/samplingFrequency)+stimDur],[-inwardORoutward*(ymax),-inwardORoutward*(ymax)],'Color',[0 0.4470 0.7410],'LineWidth',10)
end

% adding scale bar
xmaxScale = xmax;
xminScale = xmin;
line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
text(xmaxScale-(xmaxScale-xminScale)/9,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
text(xmaxScale-(xmaxScale-xminScale)/7,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/7)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
        
hold off;
movegui('north');


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

% stores sweep by sweep data
filename = strcat(fileName, '_', analysisDate, " - psc_vs_light_single - sweep_by_sweep");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(data);
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...      
    {'mouse', ...
    'date', ...
    'sweep', ...
    'discardedSweepsFromEnd', ...
    'inward(-1)ORoutward(1)', ...
    'baselineDurationInSeconds', ...
    'lightPulseAnalysisWindowInSeconds', ...
    'thresholdInDataPts', ...        
    'rsTestPulseOnsetTime', ...    
    'lightPulseDur(s)', ... 
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


% create cell array with strings for naming the mean and std of the oPSC for each light pulse
oPSCvariableNamesAmplitudeAVG = cell(1,length(lightEvokedCurrentsAmp));
oPSCvariableNamesAmplitudeSTD = cell(1,length(lightEvokedCurrentsAmp));
oPSCvariableNamesOnsetLatencyAVG = cell(1,length(allLightEvokedResponseOnsetLatencyInMilliSeconds));
oPSCvariableNamesOnsetLatencySTD = cell(1,length(allLightEvokedResponseOnsetLatencyInMilliSeconds));
oPSCvariableNamesPeakLatencyAVG = cell(1,length(allLightEvokedResponsePeakLatencyInMilliSeconds));
oPSCvariableNamesPeakLatencySTD = cell(1,length(allLightEvokedResponsePeakLatencyInMilliSeconds));
oPSCvariableNames10percentLatencyAVG = cell(1,length(allTimeTo10percentOfPeakInMilliSeconds));
oPSCvariableNames10percentLatencySTD = cell(1,length(allTimeTo10percentOfPeakInMilliSeconds));
oPSCvariableNames90percentLatencyAVG = cell(1,length(allTimeTo90percentOfPeakInMilliSeconds));
oPSCvariableNames90percentLatencySTD = cell(1,length(allTimeTo90percentOfPeakInMilliSeconds));
oPSCvariableNamesRiseTimeAVG = cell(1,length(allRiseTimeInMilliSeconds));
oPSCvariableNamesRiseTimeSTD = cell(1,length(allRiseTimeInMilliSeconds));
for lightPulse = 1:length(lightEvokedCurrentsAmp)
    oPSCvariableNamesAmplitudeAVG(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(pA)AVG')};
    oPSCvariableNamesAmplitudeSTD(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(pA)STD')};
    oPSCvariableNamesOnsetLatencyAVG(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(onsetLat_ms)AVG')};
    oPSCvariableNamesOnsetLatencySTD(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(onsetLat_ms)STD')};
    oPSCvariableNamesPeakLatencyAVG(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(peakLat_ms)AVG')};
    oPSCvariableNamesPeakLatencySTD(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(peakLat_ms)STD')};   
    oPSCvariableNames10percentLatencyAVG(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(10peakLat_ms)AVG')};
    oPSCvariableNames10percentLatencySTD(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(10peakLat_ms)STD')};
    oPSCvariableNames90percentLatencyAVG(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(90peakLat_ms)AVG')};
    oPSCvariableNames90percentLatencySTD(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(90peakLat_ms)STD')};
    oPSCvariableNamesRiseTimeAVG(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(rise_ms)AVG')};
    oPSCvariableNamesRiseTimeSTD(lightPulse) = {strcat(num2str(lightPulse), 'oPSC(rise_ms)STD')};
end

% stores cell data
filename = strcat(fileName, '_', analysisDate, " - psc_vs_light_single - cell");
fulldirectory = strcat(savefileto,'\',filename,'.xls');        
dataInCellFormat = {};
dataInCellFormat = num2cell(dataCell);
labeledData = cell2table(dataInCellFormat, 'VariableNames', ...      
    {'mouse', ...
    'date', ...
    'firstSweep', ...
    'lastSweep', ...
    'nSweeps', ... 
    'bGuess', ...
    'fastTauGuess', ...
    'slowTauGuess', ...
    'seriesResistanceAGV(Mohm)', ...
    'seriesResistanceMIN(Mohm)', ...
    'seriesResistanceMAX(Mohm)', ...   
    'baselineCurrentAVG(pA)', ...
    'baselineCurrentMIN(pA)', ...
    'baselineCurrentMAX(pA)', ...    
    oPSCvariableNamesAmplitudeAVG{:}, ...
    oPSCvariableNamesAmplitudeSTD{:}, ...
    oPSCvariableNamesOnsetLatencyAVG{:}, ...
    oPSCvariableNamesOnsetLatencySTD{:}, ...
    oPSCvariableNamesPeakLatencyAVG{:}, ...
    oPSCvariableNamesPeakLatencySTD{:}, ...
    oPSCvariableNames10percentLatencyAVG{:}, ...
    oPSCvariableNames10percentLatencySTD{:}, ...
    oPSCvariableNames90percentLatencyAVG{:}, ...
    oPSCvariableNames90percentLatencySTD{:}, ...
    oPSCvariableNamesRiseTimeAVG{:}, ...
    oPSCvariableNamesRiseTimeSTD{:}, ...
    'lightEvokedPeakCurrentAmp(pA)', ...
    'onsetLatency(ms)', ... 
    'latencyToPeak(ms)', ... 
    'latencyTo10percentOfPeak(ms)', ... 
    'latencyTo90percentOfPeak(ms', ... 
    'latencyToZeroAfterPeak(ms)', ... 
    'riseTime(ms)', ...  
    'decayFitFastTau(ms)', ...
    'decayFitSlowTau(ms)', ...
    'decayFitA', ...
    'decayFitB', ...
    'rsquare'});    
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the cell xls file')


end    