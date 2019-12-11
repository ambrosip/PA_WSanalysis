function normmonofig(obj, varargin)

[firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);

% optional arguments
numvarargs = length(varargin);
optargs = {-1 15 1.99 2 0.005 1.94 2.14 -10000 500 1 'D:\Temp\From MATLAB\test'};
optargs(1:numvarargs) = varargin;
[inwardORoutward, threshold, baselineOnsetTime, lightOnsetTime, lightDuration, xmin, xmax, ymin, ymax, rsTestPulseOnsetTime, savefileto] = optargs{:};

% variables that will be used
totalActiveChannels = 2;
sampleRate = obj.header.Acquisition.SampleRate;
lightOnsetDataPoint = lightOnsetTime*sampleRate;
lightOffDataPoint = (lightOnsetTime+lightDuration)*sampleRate;
afterLightDataPoint = (lightOnsetTime+0.01)*sampleRate;
round(afterLightDataPoint);
baselineOnsetDataPoint = baselineOnsetTime*sampleRate;
rsBaselineDataPointInterval = ((rsTestPulseOnsetTime-0.1)*sampleRate):(rsTestPulseOnsetTime*sampleRate);
rsFirstTransientDataPointInterval = (rsTestPulseOnsetTime*sampleRate):(rsTestPulseOnsetTime+0.0025)*sampleRate;
mouseNumber = getMouseNumber(obj);
experimentDate = getExperimentDate(obj);

% matrixes that will be filled
lightEvokedCurrents = [];
dataPerSweepCh1 = [];
dataPerSweepCh2 = [];
seriesResistance = [];
data = [];
allRs = [];
lightEvokedResponseLatency = [];
allLightEvokedResponseLatencyInMilliSeconds = [];

    % checking for incomplete sweeps and not analyzing incomplete sweeps (to
    % avoid this error: "Dimensions of matrices being concatenated are not
    % consistent." when calculating this: "dataPerSweepCh2 = [dataPerSweepCh2,
    % ych2]" and other concatenations.   
    
    if numel(fieldnames(obj.sweeps)) < obj.header.NSweepsPerRun
        lastSweepNumber = firstSweepNumber + numel(fieldnames(obj.sweeps)) - 2
        allSweeps = firstSweepNumber:lastSweepNumber
    end 
    
    
    % calculating all things
    for sweepNumber = allSweeps  
        [x,y] = obj.xy(sweepNumber, 1);
        
        % calculating series resistance
        rsBaselineCurrent = mean(y(rsBaselineDataPointInterval));
        rsTransientCurrent = min(y(rsFirstTransientDataPointInterval));
        dCurrent = rsTransientCurrent-rsBaselineCurrent;
        dVoltage = -5;
        seriesResistance = 1000*dVoltage/dCurrent; %mV/pA equals Gohm
        
        % normalizing data based on mean baseline current between t=1.99 (or whatever was the input time) and t=2
        y = y-mean(y(baselineOnsetDataPoint:lightOnsetDataPoint));       

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
        lightEvokedResponseLatencyInDataPoints = min(find(conv2(sign(diff(y(lightOnsetDataPoint:afterLightDataPoint))), zeros(threshold,1)+(inwardORoutward), 'valid')==threshold));
        
        % 1000 multiplication is conversion from seconds to milliseconds
        lightEvokedResponseLatencyInMilliSeconds = 1000*lightEvokedResponseLatencyInDataPoints/sampleRate;
        
        if isempty(lightEvokedResponseLatencyInMilliSeconds)
            lightEvokedResponseLatencyInMilliSeconds = NaN;
        end

        lightEvokedCurrent = min(y(lightOnsetDataPoint:afterLightDataPoint));
        
        dataPerSweepCh1 = [dataPerSweepCh1, y];
        lightEvokedCurrents = [lightEvokedCurrents, min(y(lightOnsetDataPoint:afterLightDataPoint))];  
              
        data = [data; mouseNumber, experimentDate, sweepNumber, seriesResistance, lightEvokedCurrent, lightEvokedResponseLatencyInMilliSeconds];
        allRs = [allRs, seriesResistance];
        allLightEvokedResponseLatencyInMilliSeconds = [allLightEvokedResponseLatencyInMilliSeconds, lightEvokedResponseLatencyInMilliSeconds];
        
        
    end
    
    
    
    % plotting figure with only Ch1 niceplot
        figure('name', strcat(obj.file,' (all) - niceplot - Baseline Subtracted')); % naming figure file
        hold on;
        plot(x, dataPerSweepCh1,'Color', [0.75, 0.75, 0.75]);
        plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
        rectangle('Position', [lightOnsetTime 200 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
        axis([lightOnsetTime-0.01 lightOnsetTime+0.1 ymin ymax]) % used to be -7000
        xlabel('Time (s)');
        ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
        title([obj.file ' (all) - Baseline Subtracted'],'Interpreter','none');
        set(gca,'Visible','off')
        
        xminScale = lightOnsetTime-0.01;
        xmaxScale = lightOnsetTime+0.1;
        
        % adding scale bar
        line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
        line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
        text(xmaxScale-(xmaxScale-xminScale)/10,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
        text(xmaxScale-(xmaxScale-xminScale)/8,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/7)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
                      
        hold off;
  
        
        % plotting figure with only Ch1 niceplot - AVG only
        figure('name', strcat(obj.file,' (all) - niceplot - Baseline Subtracted - AVG only')); % naming figure file
        hold on;
        plot(x, mean(dataPerSweepCh1,2),'Color','black','LineWidth',1.5); 
        rectangle('Position', [lightOnsetTime 200 lightDuration 100], 'FaceColor', [0 0.4470 0.7410], 'LineStyle', 'none')
        axis([lightOnsetTime-0.01 lightOnsetTime+0.1 ymin ymax]) % used to be -7000
        xlabel('Time (s)');
        ylabel(strcat("Baseline Subtracted ", obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorChannelName, ' (', obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits, ')'));
        title([obj.file ' (all) - Baseline Subtracted'],'Interpreter','none');
        set(gca,'Visible','off')
        
        xminScale = lightOnsetTime-0.01;
        xmaxScale = lightOnsetTime+0.1;
        
        % adding scale bar
        line([xmaxScale-(xmaxScale-xminScale)/11,xmaxScale],[ymin,ymin],'Color','k')
        line([xmaxScale,xmaxScale],[ymin,ymin+((ymax-ymin)/7)],'Color','k')
        text(xmaxScale-(xmaxScale-xminScale)/10,ymin+((ymax-ymin)/25),strcat(num2str(1000*(xmaxScale-xminScale)/11)," ms"))
        text(xmaxScale-(xmaxScale-xminScale)/8,ymin+((ymax-ymin)/10),strcat(num2str((ymax-ymin)/7)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))
                      
        hold off;
        
        
end