% EX:   s184=WSfile('m282_2019-02-12_0184-0188.h5')  load file
%       s184.plotxy(184)                             plot raw sweep 184
%       s184.plotxdy(184)                            plot deriv sweep 184
%       s184.plotxy([184,185])                       plot sweeps 184 & 185 in same figure
%       s184.plotxy([184:187])                       plot sweeps 184-187 in same figure
%       data.ff(109,0,7,10)                          prints firing frequency from sweep 109, threshold 0 mV, from 7-10 s                   

%       arrayfun(@(x) s184.plotxy(x),[184:187])      plot sweeps 184-187 in different figures

classdef WSfile
   
    properties
        header   
        sweeps
    end
    
    methods
        function obj = WSfile(fileName)
            % Need to put single quotes around file name
            % Ex: WSfile('m282_2019-02-12_0089-0090.h5')
            
            % Load file with sweeps and header using load function from WS
            loadedFile = ws.loadDataFile(fileName);
            
            % Assign header from loaded file to header property
            obj.header = loadedFile.header;
            
            % sweeps is what is left after you remove header
            obj.sweeps = rmfield(loadedFile,'header');
        end
        
        
        function loadedSweep = sweep(self, sweepNumber)
            % Ex sweepNumber: 90 for sweep_0090  
            sweepName = strcat('sweep_', num2str(sweepNumber, '%04.f')); %maybe change f (float) to integer
            loadedSweep = getfield(self.sweeps, sweepName);
            % MATLAB suggests using new notation: loadedSweep = loadedFile.(sweepName)            
        end

              
        function [x,y] = xy(self, sweepNumber) 
            % Raw y values from sweep
            % Ex sweepNumber: 90 for sweep_0090 
            samplingFrequency = self.header.Acquisition.SampleRate;
            sweepData = self.sweep(sweepNumber);
            
            % NOTE that you're getting data from Channel 1 - add new
            % variable if you want to choose channel
            rawY = sweepData.analogScans(:,1);
            
            % Need to scale y values according to Multiclamp gain
            y = rawY*self.header.Acquisition.AnalogChannelScales(1)*100;
            
            sweepDuration = size(y,1)/samplingFrequency;
            x = linspace(0,sweepDuration,size(y,1))';
        end
        
        
        function plotxy(self, sweepNumbers)
            figure;
            hold on;
            for sweepNumber = sweepNumbers
                [x,y] = self.xy(sweepNumber);
                plot(x,y)
            end
            hold off;
        end
          
        
        function [x,dy] = xdy(self, sweepNumber)
            % Derivative of y values from sweep
            samplingFrequency = self.header.Acquisition.SampleRate;
            sweepData = self.sweep(sweepNumber);
            rawY = sweepData.analogScans(:,1);
            y = rawY*self.header.Acquisition.AnalogChannelScales(1)*100;
            dy = diff(y);
            
            duration = size(y,1)/samplingFrequency;
            x = linspace(0,duration,size(y,1))';
            x(end) = [];
        end
        
        
        function plotxdy(self, sweepNumber)
            figure
            [x,y] = self.xdy(sweepNumber);
            plot(x,y)
        end
        
        
        function [pks,locs,w,p] = peaks(self, sweepNumber, MinPeakHeight)
            % MinPeakHeight is the threshold value above which function
            % will find peaks.
            
            [x,y] = self.xy(sweepNumber);
            [pks,locs,w,p] = findpeaks(y,x,'MinPeakHeight',MinPeakHeight);
            % pks   peak y value
            % locs  peak x value
            % w     peak half-width
            % p     peak amplitude
        end
        
        
        function plotp(self,sweepNumber,MinPeakHeight)
            
            [x,y] = self.xy(sweepNumber);
            [pks,locs,w,p] = findpeaks(y,x,'MinPeakHeight',MinPeakHeight);
            
            self.plotxy(sweepNumber);
            hold on;
            plot(locs,pks,'o');
            hold off;
        end
        
        
        function firingFrequency = ff(self, sweepNumber, MinPeakHeight, timeStart, timeEnd)
            [pks,locs,w,p] = self.peaks(sweepNumber, MinPeakHeight);
            
            % l = locs(locs>=timeStart && locs<= timeEnd) ;
            % numPeaks = length(l);
            
            numPeaks = sum((locs>=timeStart & locs<=timeEnd));
            firingFrequency = numPeaks/(timeEnd-timeStart);          
        end
    end
end