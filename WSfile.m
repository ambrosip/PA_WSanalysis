% EX:   s184=WSfile('m282_2019-02-12_0184-0188.h5')  load file
%       s184.plotxy(184)                             plot raw sweep 184
%       s184.plotxdy(184)                            plot deriv sweep 184
%       s184.plotxy([184,185])                       plot sweeps 184 & 185 in same figure
%       s184.plotxy([184:187])                       plot sweeps 184-187 in same figure
%       data.ff(109,0,7,10)                          prints firing frequency from sweep 109, threshold 0 mV, from 7-10 s     

%       s96.plotps(99,0,5)                           plots sweep 99, threshold 0, light onset at 5 s

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
        
        
        % varargin: variable argument input (axes range)
        function plotxy(self, sweepNumbers,varargin)
            figure;
            hold on;
            for sweepNumber = sweepNumbers
                [x,y] = self.xy(sweepNumber);
                plot(x,y)
            end
            hold off;
            
            % only want 4 optional inputs at most
            numvarargs = length(varargin);
            if numvarargs > 4
                error('myfuns:somefun2Alt:TooManyInputs', ...
                    'requires at most 4 optional inputs');
            end
            
            % set defaults for optional inputs - auto axes
            optargs = {-inf inf -inf inf};
            
            % now put these defaults into the valuesToUse cell array, 
            % and overwrite the ones specified in varargin.
            optargs(1:numvarargs) = varargin;
            
            % Place optional args in memorable variable names
            [xmin, xmax, ymin, ymax] = optargs{:};
            axis([xmin xmax ymin ymax])
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
        
        
        % plots selected peaks: 3 APs prior to light pulse, all APs during
        % light pulse, and 3 APs after light pulse. Displays firing
        % frequency as 1/ISI (ISI: Inter Spike Interval).
        function plotps(self,sweepNumber,MinPeakHeight,LightOnsetTime)
            
            [x,y] = self.xy(sweepNumber);
            [pks,locs,w,p] = findpeaks(y,x,'MinPeakHeight',MinPeakHeight);
            
            % Save locs and pks values in different variable names that 
            % will be editted differently
            locsLight = locs;
            locsBaseline1 = locs;
            locsBaseline2 = locs;
            
            pksLight = pks;
            pksBaseline1 = pks;
            pksBaseline2 = pks;
            
            % indices are integers. Need to find indices at which the value
            % of locs is between 15 and 16. Remember that locs stores the
            % time stamp at which a peak was detected
            indicesLight = find(locs<LightOnsetTime | locs>(LightOnsetTime+1));
            % deletes all peaks that happen before or after the light pulse
            locsLight(indicesLight) = [];
            pksLight(indicesLight) = [];
            
            indicesBaseline1 = find(locs>=LightOnsetTime);
            locsBaseline1(indicesBaseline1) = [];
            pksBaseline1(indicesBaseline1) = [];
            locsBaseline1(1:length(locsBaseline1)-3)=[];
            pksBaseline1(1:length(pksBaseline1)-3)=[];
            
            indicesBaseline2 = find(locs<=(LightOnsetTime+1.75));
            locsBaseline2(indicesBaseline2) = [];
            pksBaseline2(indicesBaseline2) = [];
            locsBaseline2(4:length(locsBaseline2))=[];
            pksBaseline2(4:length(pksBaseline2))=[];
            
            inverseISIforLight = 1/mean(diff(locsLight));
            inverseISIforBaseline1 = 1/mean(diff(locsBaseline1));
            inverseISIforBaseline2 = 1/mean(diff(locsBaseline2));
            
            self.plotxy(sweepNumber);
            hold on;
            plot(locsLight,pksLight,'o','color','red');
            plot(locsBaseline1,pksBaseline1,'o','color','blue');
            plot(locsBaseline2,pksBaseline2,'o','color','blue');
            rectangle('Position',[LightOnsetTime,-85,1,145],'FaceColor', [0 0 1 0.05],'LineStyle','none');
            hold off;
            axis([LightOnsetTime-2 LightOnsetTime+3 -85 45]);
            xlabel('Time (s)');
            ylabel('Voltage (mV)');
            text(LightOnsetTime-1.75,-75,'1/ISI (Hz)');
            text(LightOnsetTime+0.25,-75,num2str(inverseISIforLight),'color','red');
            text(LightOnsetTime-0.75,-75,num2str(inverseISIforBaseline1),'color','blue');
            text(LightOnsetTime+2,-75,num2str(inverseISIforBaseline2),'color','blue');
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