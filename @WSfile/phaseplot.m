function phaseplot(obj)
    
    % finding sweep numbers from file name
    if length(obj.file) == 28
        firstSweepNumber = str2num(obj.file(end-11:end-8));
        lastSweepNumber = str2num(obj.file(end-6:end-3));
        
    else
        firstSweepNumber = str2num(obj.file(end-6:end-3));
        lastSweepNumber = str2num(obj.file(end-6:end-3));
        
    end
    
    allSweeps = firstSweepNumber:lastSweepNumber;
    samplingFrequency = obj.header.Acquisition.SampleRate;
    
    % plotting individual figures for all sweeps
    for sweepNumber = allSweeps
        figure('name', strcat(obj.file,' (',num2str(sweepNumber),') Phase Plot' ));
        [x,y] = obj.xy(sweepNumber);
        [x,dy] = obj.xdy(sweepNumber);
        y(end)=[];
        plot(y,dy*samplingFrequency/1000)
        xlabel('mV');
        ylabel('mV/ms');
        title([obj.file ' (' num2str(sweepNumber) ')'],'Interpreter','none');
    end
end