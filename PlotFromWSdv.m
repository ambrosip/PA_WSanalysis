% Plot derivative

% Need to put single quotes around file name and mysweep.
% Ex: PlotFromWS('m282_2019-02-12_0089-0090.h5','sweep_0090')

function PlotWSdv = PlotFromWSdv(file,mysweep)

s = ws.loadDataFile(file);

samplingFrequency = s.header.Acquisition.SampleRate;
sweepStruct = getfield(s,mysweep);
sweepData = sweepStruct.analogScans;
sweepDuration = size(sweepData,1)/samplingFrequency;

diffData = diff(sweepData);

% Need to scale x axis accodring to sweep duration and sampling frequency
xAxis = linspace(0,sweepDuration,size(sweepData,1))';

% Need to remove last element from xAxis
xAxis(end) = [];

PlotWSdv = figure;
plot(xAxis,diffData);

end