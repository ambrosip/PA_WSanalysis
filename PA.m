classdef PA
   
    properties
        
           s
           samplingFrequency
           sweepStruct
           sweepData 
           sweepDuration 
           diffData
           
    end
    
    methods
        
        function obj = PA(file, mysweep)      
            obj.s = ws.loadDataFile(file);
            obj.samplingFrequency = obj.s.header.Acquisition.SampleRate;
            obj.sweepStruct = getfield(obj.s,mysweep);
            obj.sweepData = obj.sweepStruct.analogScans;
            obj.sweepDuration = size(obj.sweepData,1)/obj.samplingFrequency;
            obj.diffData = diff(obj.sweepData);
        end

        function PlotWS = PlotFromWS(self)
            % Plot raw data
            % Need to put single quotes around file name and mysweep.
            % Ex: PlotFromWS('m282_2019-02-12_0089-0090.h5','sweep_0090')

            % Need to scale x axis accodring to sweep duration and sampling frequency
            xAxis = linspace(0,self.sweepDuration,size(self.sweepData,1))';

            PlotWS = figure;
            plot(xAxis,self.sweepData);
        end
        
        function PlotWSdv = PlotFromWSdv(self)
            % Plot derivative
            % Need to put single quotes around file name and mysweep.
            % Ex: PlotFromWS('m282_2019-02-12_0089-0090.h5','sweep_0090')
           
            % Need to scale x axis accodring to sweep duration and sampling frequency
            xAxis = linspace(0,self.sweepDuration,size(self.sweepData,1))';

            % Need to remove last element from xAxis
            xAxis(end) = [];

            PlotWSdv = figure;
            plot(xAxis,self.diffData);

        end
        
    end
end