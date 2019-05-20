function [seriesResistance, firstSweepNumber] = calculateRs(obj)

    % finding sweep numbers from file name
    if length(obj.file) == 28
        firstSweepNumber = str2num(obj.file(end-11:end-8));
        lastSweepNumber = str2num(obj.file(end-6:end-3));
        
    else
        firstSweepNumber = str2num(obj.file(end-6:end-3));
        lastSweepNumber = str2num(obj.file(end-6:end-3));
        
    end
    
    allSweeps = firstSweepNumber:lastSweepNumber;
    
    % creates empty array to store Rs for all sweeps
    seriesResistance = [];
    
    % finds Rs for all sweeps and stores values in "seriesResistance" 
    for sweepNumber = allSweeps  
        [x,y] = obj.xy(sweepNumber, 1);
        firstTransient = min(y);
        baseline = mean(y(1:900));
        dCurrent = firstTransient-baseline;
        dVoltage = -5;
        seriesResistance = [seriesResistance 1000*dVoltage/dCurrent]; %mV/pA equals Gohm
    end
    
    % stores Rs
    seriesResistance = mean(seriesResistance);
    
    % TO DO: add timestamp info - m012.s0001.header.ClockAtRunStart
    
end