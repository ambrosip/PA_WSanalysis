%{ 
DOCUMENTATION
Created: 2022 01 26
Author: PA

This function is used to analyze the kinetics of THE FIRST post-synaptic current
(PSC) in a train, either IPSCs, EPSCs, or ChR2-mediated currents.
Originally created to study the kinetics of recordings with the internal 
CsMeSo3 with ECl = -55 mV 

This code was adapted from psc_vs_light_single

INPUTS explained:
    TO DO

    - analysisWindowAfterPeakInSeconds: When there is truly just one oIPSC, I
    analyze the decay from peak to "return to zero". Hence, this variable
    is not needed. HOWEVER, when analyzing the first oIPSC in a train, I
    noticed that the current did NOT return to zero before the start of the
    next oIPSC. Hence, I decided to crop my analysis window and look for
    the exponential fit from peak to "right before the next oIPSC". Note
    that this variable is hard coded, so it needs to be adjusted based on
    the train frequency. I chose 0.04 s (or 40 ms) as the standard when
    analyzing 20 Hz trains (which have a 50 ms inter-pulse interval).

INPUTS defaults:
    % Affects data analysis:
    lightStimCh = 2;
    discardedSweeps = [];
    discardedSweepsFromEnd = 0;
    inwardORoutward = 1;    % 1 (positive) is outward; -1 (negative) in inward
    baselineDurationInSeconds = 0.01;
    lightPulseAnalysisWindowInSeconds = 0.02;
    thresholdInDataPts = 10;
    rsTestPulseOnsetTime = 1;

    % Affects decay fit:
    bGuess = -1;            % -50    % -1000     % 1
    fastTauGuess = 0.01;    % in seconds
    slowTauGuess = 0.1;     % in seconds
    analysisWindowAfterPeakInSeconds = 0.04;

    % Affects data display:
    ymin = -1000;
    ymax = 1000;

    % Affects data saving:
    savefileto = 'D:\CORONAVIRUS DATA\From MATLAB';
    
OUTPUTS:
    TO DO

ASSUMPTIONS: 
    TO DO

TO DO:
    Fix code so that latency to 10% of peak happens between onset latency
    and peak, but NOT between pulse onset and peak!!
%}

function first_psc_kinetics(obj)
%%  USER INPUT ==================================================

% Affects data analysis:
lightStimCh = 2;
discardedSweeps = [];
discardedSweepsFromEnd = 0;
inwardORoutward = -1;    % 1 (positive) is outward; -1 (negative) in inward
baselineDurationInSeconds = 0.01;
lightPulseAnalysisWindowInSeconds = 0.02;
thresholdInDataPts = 5; %%% ALERT I CHANGED THIS FROM 10 to 5 %%%
rsTestPulseOnsetTime = 1;

% Affects decay fit:
bGuess = 0;           % -50    % -1000     % 1
fastTauGuess = 0.01;    % in seconds
slowTauGuess = 0.1;     % in seconds
analysisWindowAfterPeakInSeconds = 0.04;     % 20 Hz train has 50 ms interval between pulses. I used 40 instead of 50 to account for the latencyToPeakOnset

% Affects data display:
ymin = -250;
ymax = 250;

% Affects data saving:
savefileto = 'D:\CORONAVIRUS DATA\Out of Sync\2022 01 26 MATLAB';


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
    % Focus on getting info about the FIRST oIPSC
    lightPulseStart = find(diff(ych2>1)>0);                                      % list of light pulse onset data points
    lightPulseEnd = find(diff(ych2<1)>0);                                        % list of light pulse offset data points
    lightOnsetTime = lightPulseStart(1)/samplingFrequency;                       % onset of first light pulse in seconds
    stimDur = (lightPulseEnd(end)-lightPulseStart(end))/samplingFrequency;       % duration of each light pulse in the train (s)   
    
    xmin = lightOnsetTime - 0.01;
    xmax = lightOnsetTime + 0.05;

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

    % get info from the first oIPSC - lightPulseStart(1)
    pulseOnset = lightPulseStart(1);    
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
    
    % store all currents from a cell (only 1 column; each row is a sweep)
    lightEvokedCurrentsAllSweeps = [lightEvokedCurrentsAllSweeps; lightEvokedCurrentAmp];
    allLightEvokedResponseOnsetLatencyInMilliSecondsAllSweeps = [allLightEvokedResponseOnsetLatencyInMilliSecondsAllSweeps; lightEvokedResponseOnsetLatencyInMilliSeconds];
    allLightEvokedResponsePeakLatencyInMilliSecondsAllSweeps = [allLightEvokedResponsePeakLatencyInMilliSecondsAllSweeps; lightEvokedResponsePeakLatencyInMilliSeconds];
    allTimeTo10percentOfPeakInMilliSecondsAllSweeps = [allTimeTo10percentOfPeakInMilliSecondsAllSweeps; timeTo10percentOfPeakInMilliSeconds];
    allTimeTo90percentOfPeakInMilliSecondsAllSweeps = [allTimeTo90percentOfPeakInMilliSecondsAllSweeps; timeTo90percentOfPeakInMilliSeconds];
    allRiseTimeInMilliSecondsAllSweeps = [allRiseTimeInMilliSecondsAllSweeps; riseTimeInMilliSeconds];    
    
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

% Taking the mean of all sweeps
yBaselineSubAllMean = mean(yBaselineSubAll,2);

% ASSUMPTION ALERT
% Assuming that all sweeps have the same, single o-stim
% Update: focus on FIRST oIPSC
pulseOnset = lightPulseStart(1);
afterLightDataPoint = pulseOnset + lightPulseAnalysisWindowInDataPts;

% get amplitude and location (index) of lightEvokedCurrents
% for rise time calculation, find index of 10% and 90% of peak
% look for a threshold cross instead of an exact match - find(a<=b) instead of find(a==b) 
% outward current is +1
% inward current is 0 or -1
% for decay tau calculation, find index of return to baseline after peak
if inwardORoutward == 1 
    [lightEvokedCurrentAmp, lightEvokedCurrentLoc] = max(yBaselineSubAllMean(pulseOnset:afterLightDataPoint));
    latencyTo10percentOfPeakInDataPoints = find(yBaselineSubAllMean(pulseOnset:afterLightDataPoint) >= 0.1 * lightEvokedCurrentAmp, 1) + pulseOnset;
    latencyTo90percentOfPeakInDataPoints = find(yBaselineSubAllMean(pulseOnset:afterLightDataPoint) >= 0.9 * lightEvokedCurrentAmp, 1) + pulseOnset;
    latencyToPeakInDataPoints = lightEvokedCurrentLoc + pulseOnset;
    latencyToZeroAfterPeakInDataPoints = find(round(yBaselineSubAllMean(latencyToPeakInDataPoints:end)) <= 0, 1) + latencyToPeakInDataPoints;
else
    [lightEvokedCurrentAmp, lightEvokedCurrentLoc] = min(yBaselineSubAllMean(pulseOnset:afterLightDataPoint));
    latencyTo10percentOfPeakInDataPoints = find(yBaselineSubAllMean(pulseOnset:afterLightDataPoint) <= 0.1 * lightEvokedCurrentAmp, 1) + pulseOnset;
    latencyTo90percentOfPeakInDataPoints = find(yBaselineSubAllMean(pulseOnset:afterLightDataPoint) <= 0.9 * lightEvokedCurrentAmp, 1) + pulseOnset;
    latencyToPeakInDataPoints = lightEvokedCurrentLoc + pulseOnset;
    latencyToZeroAfterPeakInDataPoints = find(round(yBaselineSubAllMean(latencyToPeakInDataPoints:end)) >= 0, 1) + latencyToPeakInDataPoints;
end

% find onset latency
onsetLatencyInDataPoints = min(find(conv2(sign(diff(yBaselineSubAllMean(pulseOnset:afterLightDataPoint))), zeros(thresholdInDataPts,1)+(inwardORoutward), 'valid')==thresholdInDataPts));
onsetLatencyInDataPoints = onsetLatencyInDataPoints + pulseOnset;

% to stop code from breaking when it can detect onset latency in individual
% sweeps but not in the mean:
% I had to add this to deal with file m255.s0037 (recorded on 20211203) cuz
% it had 2 sweeps with good oIPSC but 1 sweep with bad oIPSC, so that the
% code found onset latency in 2/3 sweeps, but failed to find onset latency
% in the mean. So I tool the average of the onset latencies it did find and
% plotted that in the figure later.
if isempty(onsetLatencyInDataPoints)
    onsetLatencyInDataPoints = mean(allLightEvokedResponseOnsetLatencyInMilliSecondsAllSweeps, 'omitnan')*samplingFrequency/1000;
    onsetLatencyInDataPoints = onsetLatencyInDataPoints + pulseOnset;
end

% preparing to fit double-term exponential model to the decay data - create y subset
% xmaxFit used to be latencyToZeroAfterPeakInDataPoints, but in oIPSC
% trains the next oIPSC happens before the current can return to zero!
xminFit = latencyToPeakInDataPoints;
xmaxFit = latencyToPeakInDataPoints + analysisWindowAfterPeakInSeconds*samplingFrequency;
yBaselineSubAllMeanSubset = yBaselineSubAllMean;
yBaselineSubAllMeanSubset(xmaxFit:end) = [];
yBaselineSubAllMeanSubset(1:xminFit) = [];

% preparing to fit double-term exponential model to the decay data - create x subset
xSubset = x;
xSubset(xmaxFit:end) = [];
xSubset(1:xminFit) = [];
xSubset = round(xSubset,6) - round(xminFit/samplingFrequency,6);

% selecting a double exponential for decay fitting
% I was using the equation below but boyfriend said that a fitting
% algorithm runs better when you multiply instead of divide for - something
% to do with dividing by zero problems. So I substituted the Taus for their
% inverse.
% g = fittype('a*exp(-x/fastTau) + b*exp(-x/slowTau)');
g = fittype('a*exp(-x*invFastTau) + b*exp(-x*invSlowTau)');

% fit double-term exponential model to the decay data 
% If you don't give MATLAB a StartPoint, it fails horribly
% StartPoint order: [a, b, fastTau, slowTau]
% You can check the order by using this line:
% coefficientNames = coeffnames(g)
% I guessed StartPoint values by manually adjusting a double exponential
% The single exponential was NOT giving me a good fit
% I added the lower bound to stop getting negative slowTaus.
% I can also add an upper bound with the line below:
% 'Upper', [lightEvokedCurrentAmp, lightEvokedCurrentAmp, Inf, Inf]
[decayFit,gof,output] = fit(xSubset,yBaselineSubAllMeanSubset,g,...
    'StartPoint',[lightEvokedCurrentAmp, bGuess, 1/fastTauGuess, 1/slowTauGuess],...
    'Lower', [-Inf,-Inf,0,0]);

% save coefficients a and b
decayFitA = decayFit.a
decayFitB = decayFit.b

% % For troubleshooting and estimating StartPoints:
% figure;
% plot(decayFit);
% figure;
% hold on;
% plot(xSubset+round(latencyToPeakInDataPoints/samplingFrequency,6), decayFit(xSubset))
% plot(x, yBaselineSubAllMean);
% plot(xSubset+round(latencyToPeakInDataPoints/samplingFrequency,6), lightEvokedCurrentAmp*exp(-xSubset/0.05));
% hold off;

% convert data points to ms
onsetLatencyInMilliSeconds = 1000 * (onsetLatencyInDataPoints - pulseOnset) / samplingFrequency;
latencyToPeakInMilliSeconds = 1000 * (latencyToPeakInDataPoints - pulseOnset) / samplingFrequency;
latencyToZeroAfterPeakInMilliSeconds = 1000 * (latencyToZeroAfterPeakInDataPoints - pulseOnset) / samplingFrequency;
latencyTo10percentOfPeakInMilliSeconds = 1000 * (latencyTo10percentOfPeakInDataPoints - pulseOnset) / samplingFrequency;
latencyTo90percentOfPeakInMilliSeconds = 1000 * (latencyTo90percentOfPeakInDataPoints - pulseOnset) / samplingFrequency;

% convert s to ms
% decayFitFastTau = decayFit.fastTau * 1000
% decayFitSlowTau = decayFit.slowTau * 1000
decayFitFastTau = (1/decayFit.invFastTau) * 1000
decayFitSlowTau = (1/decayFit.invSlowTau) * 1000

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
                      

%% PLOT - rs =====================================================

figure('name', strcat(fileName, " ", analysisDate, ' - first_psc_kinetics - Rs all')); % naming figure file
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
figure('name', strcat(fileName, " ", analysisDate, ' - first_psc_kinetics - Niceplot all colors'));
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
title([fileName ' - first_psc_kinetics - niceplot all colors'],'Interpreter','none');
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
line([xmaxScale-(xmaxScale-xminScale)/10,xmaxScale],[ymin,ymin],'Color','k')
line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/5)],'Color','k')
text(xmaxScale-(xmaxScale-xminScale)/8,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/10)," ms"))
text(xmaxScale-(xmaxScale-xminScale)/6,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/5)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
hold off;
movegui('southeast');


%% PLOT - niceplot of all sweeps and/or MEAN - BLACK ================================

% plotting niceplot     
figure('name', strcat(fileName, " ", analysisDate, ' - first_psc_kinetics - niceplot all'));
hold on;
plot(x, yBaselineSubAll,'Color',[0, 0, 0, 0.25]);
plot(x, mean(yBaselineSubAll,2),'Color','black','LineWidth',0.7); 
line([xmin xmax],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
axis([xmin xmax ymin ymax]);
xlabel('Time (s)');
ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
title([fileName ' - first_psc_kinetics - niceplot all'],'Interpreter','none');
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
line([xmaxScale-(xmaxScale-xminScale)/10,xmaxScale],[ymin,ymin],'Color','k')
line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/5)],'Color','k')
text(xmaxScale-(xmaxScale-xminScale)/9,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/10)," ms"))
text(xmaxScale-(xmaxScale-xminScale)/7,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/5)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
hold off;
movegui('south');


%% PLOT - subtracted baseline current =====================================================

figure('name', strcat(fileName, " ", analysisDate, ' - first_psc_kinetics - baseline current')); % naming figure file
plot(allSweeps, baselineCurrentAll,'-o');
axis([allSweeps(1) inf -300 300])
ylabel('Subtracted Baseline Current (pA)');
xlabel('Sweeps');
title([obj.file ' rs'],'Interpreter','none');
movegui('southwest');


%% PLOT - niceplot with kinetics - BLACK ================================

% plotting niceplot     
figure('name', strcat(fileName, " ", analysisDate, ' - first_psc_kinetics - Niceplot kinetics black'));
hold on;
% plot(x, yBaselineSubAll,'Color',[0, 0, 0, 0.25]);
plot(x, mean(yBaselineSubAll,2),'Color','black','LineWidth',0.7); 

%%%% ALERT ALERT on 2022 03 03, I added the "-1.5" to align the lined to
%%%% actual points in the data. I'm not sure I did the right thing. I used
%%%% the peak data as my guide.
plot(xSubset+round((latencyToPeakInDataPoints-1.5)/samplingFrequency,6), decayFit(xSubset))
plot((latencyToPeakInDataPoints-1.5)/samplingFrequency, lightEvokedCurrentAmp, 'o', 'Color', 'red');
line([(onsetLatencyInDataPoints-1.5)/samplingFrequency (onsetLatencyInDataPoints-1.5)/samplingFrequency], [ymin ymax], 'Color',[0.5 0.5 0.5],'LineStyle','--');
line([(latencyTo10percentOfPeakInDataPoints-1.5)/samplingFrequency (latencyTo10percentOfPeakInDataPoints-1.5)/samplingFrequency], [ymin ymax], 'Color',[0.5 0.5 0.5],'LineStyle','--');
line([(latencyTo90percentOfPeakInDataPoints-1.5)/samplingFrequency (latencyTo90percentOfPeakInDataPoints-1.5)/samplingFrequency], [ymin ymax], 'Color',[0.5 0.5 0.5],'LineStyle','--');
line([latencyToZeroAfterPeakInDataPoints latencyToZeroAfterPeakInDataPoints], [ymin ymax], 'Color',[0.5 0.5 0.5],'LineStyle','--');
line([xmin xmax],[0, 0],'Color',[0.5 0.5 0.5],'LineStyle','--')
axis([xmin xmax ymin ymax]);
xlabel('Time (s)');
ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
title([fileName ' - first_psc_kinetics - niceplot kinetics'],'Interpreter','none');
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
line([xmaxScale-(xmaxScale-xminScale)/10,xmaxScale],[ymin,ymin],'Color','k')
line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/5)],'Color','k')
text(xmaxScale-(xmaxScale-xminScale)/9,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/10)," ms"))
text(xmaxScale-(xmaxScale-xminScale)/7,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/5)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
        
hold off;
movegui('north');


%% EXPORTING XLS files ==========================================

% stores sweep by sweep data
filename = strcat(fileName, '_', analysisDate, " - first_psc_kinetics - sweep_by_sweep");
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
    'oPSC(pA)', ...
    'oPSC(onsetLat_ms)', ...
    'oPSC(peakLat_ms)', ...
    'oPSC(10peakLat_ms)', ... 
    'oPSC(90peakLat_ms)', ...
    'oPSC(riseTime_ms)'});    
writetable(labeledData, fulldirectory, 'WriteMode', 'overwritesheet');
disp('I saved the sweep_by_sweep xls file')


% stores cell data
filename = strcat(fileName, '_', analysisDate, " - first_psc_kinetics - cell");
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
    'oPSC(pA)AVG', ...
    'oPSC(pA)STD', ...
    'oPSC(onsetLat_ms)AVG', ...
    'oPSC(onsetLat_ms)STD', ...
    'oPSC(peakLat_ms)AVG', ...
    'oPSC(peakLat_ms)STD', ...
    'oPSC(10peakLat_ms)AVG', ...
    'oPSC(10peakLat_ms)STD', ...
    'oPSC(90peakLat_ms)AVG', ...
    'oPSC(90peakLat_ms)STD', ...
    'oPSC(rise_ms)AVG', ...
    'oPSC(rise_ms)STD', ...
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