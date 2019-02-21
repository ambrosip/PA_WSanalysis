% Plot raw data

% Need to put single quotes around file name and mysweep.
% Ex: PlotFromWS('m282_2019-02-12_0089-0090.h5','sweep_0090')

function PlotWS = PlotFromWS(file,mysweep)

s = ws.loadDataFile(file);

samplingFrequency = s.header.Acquisition.SampleRate;
sweepStruct = getfield(s,mysweep);
sweepData = sweepStruct.analogScans;
sweepDuration = size(sweepData,1)/samplingFrequency;

% Need to scale x axis accodring to sweep duration and sampling frequency
xAxis = linspace(0,sweepDuration,size(sweepData,1))';

PlotWS = figure;
plot(xAxis,sweepData);

end

