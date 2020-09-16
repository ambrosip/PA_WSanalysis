function [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj)

 % finding sweep numbers from file name
    % if files are named in my naming convention
    if length(obj.file) == 28   
        firstSweepNumber = str2num(obj.file(end-11:end-8));
        lastSweepNumber = str2num(obj.file(end-6:end-3));
        
    % if files have one single sweep    
    else
        firstSweepNumber = str2num(obj.file(end-6:end-3));
        lastSweepNumber = str2num(obj.file(end-6:end-3));
        
    end
    
    allSweeps = firstSweepNumber:lastSweepNumber;
    
end