classdef newWSfile
    
    % Use with files recorded with new version of wavesurfer
    % By Mar 2020 only Talia's rig has this version
    
    properties
        header   
        sweeps
        file
    end
    
    methods
        
        newWSplot(obj, varargin);
        newWSnormmono(obj, varargin);
        
        function obj = newWSfile(fileName)                        
            % Load file with sweeps and header using load function from WS
            loadedFile = ws.loadDataFile(fileName);
            % Need to put single quotes around file name
            % Ex: WSfile('m282_2019-02-12_0089-0090.h5')
            
            % Assign header from loaded file to header property
            obj.header = loadedFile.header;
            
            % sweeps is what is left after you remove header
            obj.sweeps = rmfield(loadedFile,'header');
            
            obj.file = fileName;            
        end        

        
        function loadedSweep = sweep(self, sweepNumber)
            % Ex sweepNumber: 90 for sweep_0090  
            sweepName = strcat('sweep_', num2str(sweepNumber, '%04.f')); %maybe change f (float) to integer
            loadedSweep = getfield(self.sweeps, sweepName);
            % MATLAB suggests using new notation: loadedSweep = loadedFile.(sweepName)            
        end

        
        function [x,y] = xy(self, sweepNumber, channel)             
            % optional arguments: channel
            if nargin < 3
                selectedChannel = 1;
            else
                selectedChannel = channel;
            end
            % Raw y values from sweep
            % Ex sweepNumber: 90 for sweep_0090 
            samplingFrequency = self.header.AcquisitionSampleRate; % This line changed from original WS code
            sweepData = self.sweep(sweepNumber);        
            rawY = sweepData.analogScans(:,selectedChannel);
            y = rawY;         
            sweepDuration = size(y,1)/samplingFrequency;
            x = linspace(0,sweepDuration,size(y,1))';
        end

        
        function experimentDate = getExperimentDate(obj)
            experimentDate = str2num(strcat(obj.file(6:9),obj.file(11:12),obj.file(14:15)));
        end
 
        
        function mouseNumber = getMouseNumber(obj)
            mouseNumber = str2num(obj.file(2:4));
        end
 
        
        function [firstSweepNumber, lastSweepNumber, allSweeps] = getSweepNumbers(obj)
            % finding sweep numbers from file name
            if length(obj.file) == 28
                firstSweepNumber = str2num(obj.file(end-11:end-8));
                lastSweepNumber = str2num(obj.file(end-6:end-3));
            else
                firstSweepNumber = str2num(obj.file(end-6:end-3));
                lastSweepNumber = str2num(obj.file(end-6:end-3));
            end
            allSweeps = firstSweepNumber:lastSweepNumber;           
        end          
    end
end
