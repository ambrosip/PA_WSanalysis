function filteredplot(obj, varargin)

    %% BEFORE USING THIS FUNCTION
    
    % Add the folder where your h5 files to MATLAB's path
    
    % Load the h5 files using the fuction WSfileFromDir. Name the object
    % using the mouse number and make sure there are 3 numbers. If the
    % mouse number is m23, name it as m023.   
    % EX: m155=WSfileFromDir('E:\Priscilla\Ephys\20190808 dms chr2 habit')
    % TIP: don't load the files directly from the server because that takes
    % a million years. 
    
    % Choose which files you want to analyze - you will have to tell matlab
    % what is the first sweep number of the file you need analyzed.
    
    % Run the function using the optional arguments (or not) - feel free to 
    % edit the default values of the optional arguments to speed up your analysis
    % EX: firingrate(m155.s0201)
    % EX: firingrate(m155.s0201,10,1,5) - this one sets an optional
    % threshold and start-end points for firing rate analysis

    %% BEHOLD THE CODE
    
    % optional arguments
    % set defaults for optional inputs 
    optargs = {1 0.2 1000 5 25 -800 -800 100 1 'D:\Temp\From MATLAB 2020'}; % used for minis data
%     optargs = {1 2 1000 15 20 -400 -200 100 1 'D:\Temp\From MATLAB 2020'};
%     optargs = {1 2 1000 0 30 -300 -200 100 1 'D:\Temp\From MATLAB 2020'};
    
    % overwrite defaults with values spclose ecified in varargin
    numvarargs = length(varargin);
    optargs(1:numvarargs) = varargin;
    
    % place optional args in memorable variable names
    [filterORnot, highpassThreshold, lowpassThreshold, xmin, xmax, yminRaw, ymin, ymax, rsTestPulseOnsetTime, savefileto] = optargs{:};

    samplingFrequency = obj.header.Acquisition.SampleRate;
    sampleRate = samplingFrequency;
    rsBaselineDataPointInterval = ((rsTestPulseOnsetTime-0.1)*sampleRate):(rsTestPulseOnsetTime*sampleRate);
    rsFirstTransientDataPointInterval = (rsTestPulseOnsetTime*sampleRate):(rsTestPulseOnsetTime+0.0025)*sampleRate;
    allRs = [];

    [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj);

    % checking for incomplete sweeps and not analyzing them
    if numel(fieldnames(obj.sweeps)) < obj.header.NSweepsPerRun
        lastSweepNumber = firstSweepNumber + numel(fieldnames(obj.sweeps)) - 2;
        allSweeps = firstSweepNumber:lastSweepNumber;
    end 
    
    % plotting one figure per sweep that shows whole sweep
    for sweepNumber = allSweeps
        figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - raw')); % naming figure file
        [x,y] = obj.xy(sweepNumber, 1);
        plot(x,y,'k','LineWidth',1);
        axis([xmin xmax yminRaw ymax]);           
        xlabel('Time (s)');
        ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
        title([' (' num2str(sweepNumber) ') - raw'],'Interpreter','none');
        set(gca,'Visible','off')        
        % adding scale bar
        line([xmin,xmin+(xmax-xmin)/10],[yminRaw,yminRaw],'Color','k')
        line([xmin,xmin],[yminRaw,yminRaw+((ymax-yminRaw)/10)],'Color','k')
        text(xmin+(xmax-xmin)/130,yminRaw+((ymax-yminRaw)/30),strcat(num2str((xmax-xmin)/10)," s"))
        text(xmin+(xmax-xmin)/130,yminRaw+((ymax-yminRaw)/12),strcat(num2str((ymax-yminRaw)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))   
        % adding line at 0 pA
        line([xmin xmax],[0,0],'Color','black','LineStyle','--')
        
        % calculating series resistance
        rsBaselineCurrent = mean(y(rsBaselineDataPointInterval));
        rsTransientCurrent = min(y(rsFirstTransientDataPointInterval));
        dCurrent = rsTransientCurrent-rsBaselineCurrent;
        dVoltage = -5;
        seriesResistance = 1000*dVoltage/dCurrent; %mV/pA equals Gohm
        allRs = [allRs, seriesResistance];
        
        if filterORnot
           figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - low and highpass filtered')); % naming figure file
           yFilteredOnce = highpass(y, highpassThreshold,samplingFrequency);
           yFilteredTwice = lowpass(yFilteredOnce,lowpassThreshold,samplingFrequency);
           plot(x,yFilteredTwice,'k','LineWidth',1);
           axis([xmin xmax ymin ymax]);           
           xlabel('Time (s)');
           ylabel(obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits);
           title([' (' num2str(sweepNumber) ') - low and highpass filtered'],'Interpreter','none');
           set(gca,'Visible','off')        
           % adding scale bar
           line([xmin,xmin+(xmax-xmin)/10],[ymin,ymin],'Color','k')
           line([xmin,xmin],[ymin,ymin+((ymax-ymin)/10)],'Color','k')
           text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/30),strcat(num2str((xmax-xmin)/10)," s"))
           text(xmin+(xmax-xmin)/130,ymin+((ymax-ymin)/12),strcat(num2str((ymax-ymin)/10)," ",obj.header.Ephys.ElectrodeManager.Electrodes.element1.MonitorUnits))   
        end        
    end
    
    % plot rs
    figure('name', strcat(obj.file,' (all) - rs')); % naming figure file
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
        
    saveAllFigs(savefileto);
    disp('I saved the figures!')

end