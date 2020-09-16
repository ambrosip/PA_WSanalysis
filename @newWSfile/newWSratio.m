% Code adapted from newWSppr
% Use this code for AMPA to NMDA ratio analysis
% This code takes two objects as input:
% objTotal   total current (before D-AP5 washon)
% objAmpa    AMPA current (after D-AP5 washon)
% ex: newWSratio(m051.s0080,m051.s0089)

function newWSratio(objTotal, objAmpa, varargin)

[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(objTotal);
[firstSweepNumber1, lastSweepNumber1, allSweeps1] = getSweepNumbers(objAmpa);

% optional arguments
numvarargs = length(varargin);
optargs = {1 0.015 10 0.2 0.2 0.9 -50 300 0.1 'D:\CORONAVIRUS DATA\From MATLAB'}; 
optargs(1:numvarargs) = varargin;
[inwardORoutward, analysisTimeWindow, threshold, baselineOnsetTime, xmin, xmax, ymin, ymax, rsTestPulseOnsetTime, savefileto] = optargs{:};

% variables that will be used
totalActiveChannels = 2;
sampleRate = objTotal.header.AcquisitionSampleRate;
baselineOnsetDataPoint = baselineOnsetTime*sampleRate;
rsBaselineDataPointInterval = ((rsTestPulseOnsetTime-0.05)*sampleRate):(rsTestPulseOnsetTime*sampleRate);
rsFirstTransientDataPointInterval = (rsTestPulseOnsetTime*sampleRate):(rsTestPulseOnsetTime+0.0025)*sampleRate;
mouseNumber = getMouseNumber(objTotal);
experimentDate = getExperimentDate(objTotal);

% matrixes that will be filled
allPeakCurrentTotal = [];
allPeakCurrentAmpa = [];
dataPerSweepCh1Total = [];
dataPerSweepCh2Total = [];
dataPerSweepCh1Ampa = [];
dataPerSweepCh2Ampa = [];
dataTotal = [];
dataAmpa = [];
allRsTotal = [];
allRsAmpa = [];
allOnsetLatencyInMilliSecondsTotal = [];
allPeakLatencyInMilliSecondsTotal = []; 
allOnsetLatencyInMilliSecondsAmpa = [];
allPeakLatencyInMilliSecondsAmpa = []; 

% checking for incomplete sweeps for objTotal
% and not analyzing incomplete sweeps (to
% avoid this error: "Dimensions of matrices being concatenated are not
% consistent." when calculating this: "dataPerSweepCh2 = [dataPerSweepCh2,
% ych2]" and other concatenations.       
if numel(fieldnames(objTotal.sweeps)) < objTotal.header.NSweepsPerRun
    lastSweepNumber = firstSweepNumber + numel(fieldnames(objTotal.sweeps)) - 2;
    allSweeps = firstSweepNumber:lastSweepNumber;
end

% checking for incompete sweeps for objAmpa
if numel(fieldnames(objAmpa.sweeps)) < objAmpa.header.NSweepsPerRun
    lastSweepNumber1 = firstSweepNumber1 + numel(fieldnames(objAmpa.sweeps)) - 2;
    allSweeps1 = firstSweepNumber1:lastSweepNumber1;
end
    
%% getting all the data for objTotal
for sweepNumber = allSweeps

    % getting light stim data from ch2
    [xch2,ych2Total] = objTotal.xy(sweepNumber, 2);     
    dataPerSweepCh2Total = [dataPerSweepCh2Total, ych2Total];

    % finding info about light stim        
    lightPulseStartList = find(diff(ych2Total>1)>0);
    lightPulseEndList = find(diff(ych2Total<1)>0);
    lightPulseOnsetDataPoint = lightPulseStartList(1);
    lightPulseOnsetTime = lightPulseStartList(1)/sampleRate;                           % in seconds
    stimDur = (lightPulseEndList(end)-lightPulseStartList(end))/sampleRate;            % duration of the light pulse (s)
    afterLightDataPoint = (lightPulseOnsetTime+analysisTimeWindow)*sampleRate;

    [xTotal,yTotal] = objTotal.xy(sweepNumber, 1);

    % calculating series resistance
    rsBaselineCurrentTotal = mean(yTotal(rsBaselineDataPointInterval));
    rsTransientCurrentTotal = min(yTotal(rsFirstTransientDataPointInterval));
    dCurrentTotal = rsTransientCurrentTotal-rsBaselineCurrentTotal;
    dVoltage = -5;
    seriesResistanceTotal = 1000*dVoltage/dCurrentTotal; %mV/pA equals Gohm      

    % normalizing data based on mean current 10 ms before baselineOnsetDataPoint
    yTotal = yTotal-mean(yTotal(baselineOnsetDataPoint-10:baselineOnsetDataPoint));

    % calculating latency 
    % rationale for finding light evoked response latency: I am looking
    % for a monotonic change of the signal for at least "x" data points 
    % (x=threshold), and I am selecting the first occurence of this monotonic change.

    % the expected direction of the monotonic change is determined
    % by the optarg "inwardORoutward". If the value is 1, the function
    % will look for a monotonic increase (outward current), and if the
    % value is -1, the function will look for a monotonic decrease
    % (inward current).

    % y(lightOnsetDataPoint:lightOffDataPoint) is the signal during light pulse

    % zeros(threshold,1)+(inwardORoutward) is a vector with "threshold" 
    % number of rows containing -1 or 1, depending on the value of the
    % optarg "inwardORoutward"

    % function diff calculates difference between data points

    % function sign only keeps info about decrease (-1) or increase (1) in signal

    % function conv2 performs convolution, looking for a sequence of
    % "threshold" points of value equal to the value stored in
    % "inwardORoutward" (-1 or 1). I think conv2 collapses all found 
    % elements into a single element of value "threshold" and moves 
    % forward in its search. So if threshold=20, inwardORoutward=-1,
    % and there are 40 elements of value -1 in sequence (one after the
    % other), conv2 will spit out a vector with 2 elements, both of
    % value 20.

    % function find looks for the index of the elements equal to threshold in the output of the convolution

    % function min looks for the min index (aka the first data point
    % that marks the start of a monotonic change in signal for "x"
    % points (x=threshold).
    onsetLatencyInDataPointsTotal = min(find(conv2(sign(diff(yTotal(lightPulseOnsetDataPoint:afterLightDataPoint))), zeros(threshold,1)+(inwardORoutward), 'valid')==threshold));
    
    % 1000 multiplication is conversion from seconds to milliseconds
    onsetLatencyInMilliSecondsTotal = 1000*onsetLatencyInDataPointsTotal/sampleRate;
    
    % to avoid errors due to empty matrices
    if isempty(onsetLatencyInMilliSecondsTotal)
        onsetLatencyInMilliSecondsTotal = NaN;
    end

    % storing normalized data
    dataPerSweepCh1Total = [dataPerSweepCh1Total, yTotal];

    % storing relevant info
    % need to check if I'm looking for an inward or outward current
    if inwardORoutward == 1
        peakCurrentTotal = max(yTotal(lightPulseOnsetDataPoint:afterLightDataPoint));
    else
        peakCurrentTotal = min(yTotal(lightPulseOnsetDataPoint:afterLightDataPoint));
    end
    peakLatencyDataPointTotal = find(yTotal(lightPulseOnsetDataPoint:afterLightDataPoint)==peakCurrentTotal);
    peakLatencyInMilliSecondsTotal = 1000*peakLatencyDataPointTotal(1,1)/sampleRate;                

    % not sure I need this?
    allPeakCurrentTotal = [allPeakCurrentTotal, peakCurrentTotal];  

    % create table that will be saved as csv file      
    dataTotal = [dataTotal;...
        mouseNumber, experimentDate, sweepNumber, seriesResistanceTotal,...
        stimDur, peakCurrentTotal,...
        onsetLatencyInMilliSecondsTotal, peakLatencyInMilliSecondsTotal];

    % create lists that will be used to plot Rs and latency x sweep
    allRsTotal = [allRsTotal, seriesResistanceTotal];
    allOnsetLatencyInMilliSecondsTotal = [allOnsetLatencyInMilliSecondsTotal, onsetLatencyInMilliSecondsTotal];
    allPeakLatencyInMilliSecondsTotal = [allPeakLatencyInMilliSecondsTotal, peakLatencyInMilliSecondsTotal];  
   
end

%% getting all the data for objAmpa
for sweepNumber1 = allSweeps1

    % getting light stim data from ch2
    [xch2,ych2Ampa] = objAmpa.xy(sweepNumber1, 2);     
    dataPerSweepCh2Ampa = [dataPerSweepCh2Ampa, ych2Ampa];

    % finding info about light stim        
    lightPulseStartList = find(diff(ych2Ampa>1)>0);
    lightPulseEndList = find(diff(ych2Ampa<1)>0);
    lightPulseOnsetDataPoint = lightPulseStartList(1);
    lightPulseOnsetTime = lightPulseStartList(1)/sampleRate;                           
    stimDur = (lightPulseEndList(end)-lightPulseStartList(end))/sampleRate;            
    afterLightDataPoint = (lightPulseOnsetTime+analysisTimeWindow)*sampleRate;

    [xAmpa,yAmpa] = objAmpa.xy(sweepNumber1, 1);

    % calculating series resistance
    rsBaselineCurrentAmpa = mean(yAmpa(rsBaselineDataPointInterval));
    rsTransientCurrentAmpa = min(yAmpa(rsFirstTransientDataPointInterval));
    dCurrentAmpa = rsTransientCurrentAmpa-rsBaselineCurrentAmpa;
    dVoltage = -5;
    seriesResistanceAmpa = 1000*dVoltage/dCurrentAmpa; %mV/pA equals Gohm      

    % normalizing data based on mean current 10 ms before baselineOnsetDataPoint
    yAmpa = yAmpa-mean(yAmpa(baselineOnsetDataPoint-10:baselineOnsetDataPoint));

    % calculating latency 
    onsetLatencyInDataPointsAmpa = min(find(conv2(sign(diff(yAmpa(lightPulseOnsetDataPoint:afterLightDataPoint))), zeros(threshold,1)+(inwardORoutward), 'valid')==threshold));
    
    % 1000 multiplication is conversion from seconds to milliseconds
    onsetLatencyInMilliSecondsAmpa = 1000*onsetLatencyInDataPointsAmpa/sampleRate;
    
    % to avoid errors due to empty matrices
    if isempty(onsetLatencyInMilliSecondsAmpa)
        onsetLatencyInMilliSecondsAmpa = NaN;
    end

    % storing normalized data
    dataPerSweepCh1Ampa = [dataPerSweepCh1Ampa, yAmpa];

    % storing relevant info
    if inwardORoutward == 1
        peakCurrentAmpa = max(yAmpa(lightPulseOnsetDataPoint:afterLightDataPoint));
    else
        peakCurrentAmpa = min(yAmpa(lightPulseOnsetDataPoint:afterLightDataPoint));
    end
    peakLatencyDataPointAmpa = find(yAmpa(lightPulseOnsetDataPoint:afterLightDataPoint)==peakCurrentAmpa);
    peakLatencyInMilliSecondsAmpa = 1000*peakLatencyDataPointAmpa(1,1)/sampleRate;                

    % not sure I need this?
    allPeakCurrentAmpa = [allPeakCurrentAmpa, peakCurrentAmpa];  

    % create table that will be saved as csv file      
    dataAmpa = [dataAmpa;...
        mouseNumber, experimentDate, sweepNumber1, seriesResistanceAmpa,...
        stimDur, peakCurrentAmpa,...
        onsetLatencyInMilliSecondsAmpa, peakLatencyInMilliSecondsAmpa];

    % create lists that will be used to plot Rs and latency x sweep
    allRsAmpa = [allRsAmpa, seriesResistanceAmpa];
    allOnsetLatencyInMilliSecondsAmpa = [allOnsetLatencyInMilliSecondsAmpa, onsetLatencyInMilliSecondsAmpa];
    allPeakLatencyInMilliSecondsAmpa = [allPeakLatencyInMilliSecondsAmpa, peakLatencyInMilliSecondsAmpa];  
   
end

%% plotting

% dataPerSweepCh1 has a column for each sweep. mean(dataPerSweepCh1,2)
% calculates the average column. mean(dataPerSweepCh1) calculates the
% wrong thing: the average row.   

% plotting niceplot showing NMDA (Total - AMPA) and AMPA currents
    % naming figure file
    figure('name', strcat(objTotal.file,' (all) - niceplot - Nmda and Ampa')); 
    hold on;
%     rawNmda = plot(xTotal, dataPerSweepCh1Total-dataPerSweepCh1Ampa,'Color', [0 0 0 0.05]);
%     rawAmpa = plot(xTotal, dataPerSweepCh1Ampa,'Color', [213/255 94/255 0/255 0.05]);
    plot(xTotal, mean(dataPerSweepCh1Total,2)-mean(dataPerSweepCh1Ampa,2),'Color',[0 0 0],'LineWidth',1.5); 
    plot(xTotal, mean(dataPerSweepCh1Ampa,2),'Color',[213/255 94/255 0/255],'LineWidth',1.5);
    rectangle('Position', [lightPulseOnsetTime ymax-20 stimDur 20], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
    axis([xmin xmax ymin ymax]);
    xlabel('Time (s)');
    ylabel(strcat("Baseline subtracted ", objTotal.header.AIChannelNames(1), ' (', objTotal.header.AIChannelUnits(1), ')'));
    set(gca,'Visible','off')
    % adding scale bar
    line([xmax-(xmax-xmin)/10,xmax],[ymax,ymax],'Color','k')
    line([xmax,xmax],[ymax,ymax-((ymax-ymin)/10)],'Color','k')
    text(xmax-(xmax-xmin)/10,ymax-((ymax-ymin)/25),strcat(num2str(1000*(xmax-xmin)/10)," ms"))
    text(xmax-(xmax-xmin)/8,ymax-((ymax-ymin)/10),strcat(num2str((ymax-ymin)/10)," ",objTotal.header.AIChannelUnits(1)))     
    hold off;

% plotting niceplot
    % naming figure file
    figure('name', strcat(objTotal.file,' (all) - niceplot - Total and Ampa')); 
    hold on;
    rawTotal = plot(xTotal, dataPerSweepCh1Total,'Color', [0 0 0 0.05]);
    rawAmpa = plot(xTotal, dataPerSweepCh1Ampa,'Color', [213/255 94/255 0/255 0.05]);
    plot(xTotal, mean(dataPerSweepCh1Total,2),'Color',[0 0 0],'LineWidth',1.5); 
    plot(xTotal, mean(dataPerSweepCh1Ampa,2),'Color',[213/255 94/255 0/255],'LineWidth',1.5);
    rectangle('Position', [lightPulseOnsetTime ymax-20 stimDur 20], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
    axis([xmin xmax ymin ymax]);
    xlabel('Time (s)');
    ylabel(strcat("Baseline subtracted ", objTotal.header.AIChannelNames(1), ' (', objTotal.header.AIChannelUnits(1), ')'));
    set(gca,'Visible','off')
    % adding scale bar
    line([xmax-(xmax-xmin)/10,xmax],[ymax,ymax],'Color','k')
    line([xmax,xmax],[ymax,ymax-((ymax-ymin)/10)],'Color','k')
    text(xmax-(xmax-xmin)/10,ymax-((ymax-ymin)/25),strcat(num2str(1000*(xmax-xmin)/10)," ms"))
    text(xmax-(xmax-xmin)/8,ymax-((ymax-ymin)/10),strcat(num2str((ymax-ymin)/10)," ",objTotal.header.AIChannelUnits(1)))     
    hold off;
    
% plot data and found peaks
    % naming figure file
    figure('name', strcat(objTotal.file,' - Baseline Subtracted - Total and Ampa')); 
    hold on;   
    % plotting currents
    plot(xTotal, dataPerSweepCh1Total(:,:),'color',[0 0 0 0.2]);
    plot(xTotal, dataPerSweepCh1Ampa(:,:),'color',[213/255 94/255 0/255 0.2]); 
    plot(lightPulseOnsetTime + dataTotal(:,8)/1000, dataTotal(:,6),'o','color','red');
    plot(lightPulseOnsetTime + dataAmpa(:,8)/1000, dataAmpa(:,6),'o','color','red');
    % plotting light stim
    rectangle('Position', [lightPulseOnsetTime ymax-20 stimDur 20], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none');  
    axis([xmin xmax ymin ymax]);
    xlabel('Time (s)');
    ylabel(strcat("Baseline subtracted ", objTotal.header.AIChannelNames(1), ' (', objTotal.header.AIChannelUnits(1), ')'));
    title([objTotal.file ' Total and AMPA'],'Interpreter','none');
    hold off;
    movegui('northwest');
    
% plot rs for objTotal and objAmpa 
    figure('name', strcat(objTotal.file,' (all) - rs - Total and Ampa')); % naming figure file
    hold on;
    plot(allSweeps, allRsTotal,'-o');
    plot(allSweeps1, allRsAmpa,'-o','Color','red');
    % plot lines marking 30% increase and 30% decrese in Rs compared to first
    % test pulse
    line([allSweeps(1) allSweeps1(end)],[allRsTotal(1)*0.7, allRsTotal(1)*0.7],'Color','black','LineStyle','--')
    line([allSweeps(1) allSweeps1(end)],[allRsTotal(1)*1.3, allRsTotal(1)*1.3],'Color','black','LineStyle','--')
    axis([allSweeps(1) inf 0 60])
    hold off;
    ylabel('Rs (M\Omega)');
    xlabel('Sweeps');
    title([objTotal.file ' Total and AMPA rs'],'Interpreter','none');
    movegui('northeast');

% testing if the same number of sweeps were done in each file (aka the 
% data matrices will have the same size) to avoid this error: "Dimensions 
% of arrays being concatenated are not consistent"
    if size(dataTotal,1) > size(dataAmpa,1)
        dataAmpa(size(dataTotal,1),:)=NaN;
    elseif size(dataTotal,1) < size(dataAmpa,1)
        dataTotal(size(dataAmpa,1),:)=NaN;
    else
        disp('oof dataAmpa and dataTotal are the same size!')
    end            
    
% save xls file with all the data
    filename = strcat(objTotal.file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - AMPA NMDA ratio",'(',num2str(inwardORoutward),')');
    fulldirectory = strcat(savefileto,'\',filename,'.xls');        
    dataInCellFormat = {};
    allData = [dataTotal,dataAmpa];
    dataInCellFormat = num2cell(allData);
    labeledData = cell2table(dataInCellFormat,'VariableNames',...
        {'mouse', 'date', 'sweep', 'seriesResistance',...
        'stimDurInSeconds', 'peakCurrentTotal',...
        'onsetLatencyInMilliSeconds', 'peakLatencyInMilliSeconds',...
        'mouseAmpa', 'dateAmpa', 'sweepAmpa', 'seriesResistanceAmpa',...
        'stimDurInSecondsAmpa', 'peakCurrentAmpa',...
        'onsetLatencyInMilliSecondsAmpa', 'peakLatencyInMilliSecondsAmpa'
        });
    writetable(labeledData,fulldirectory);
    disp('I saved it')
            
end