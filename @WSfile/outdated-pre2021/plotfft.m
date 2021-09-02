function plotfft(obj,sweepNumber,highpassThreshold,lowpassThreshold)

sweepDuration = obj.header.Acquisition.Duration;
samplingFrequency = obj.header.Acquisition.SampleRate;
lengthOfSignal = sweepDuration*samplingFrequency;

[x,y] = obj.xy(sweepNumber, 1);
Y = fft(y);
P2 = abs(Y/lengthOfSignal);
P1 = P2(1:lengthOfSignal/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = samplingFrequency*(0:(lengthOfSignal/2))/lengthOfSignal;

% butter(20,

% % % figure('name', strcat(obj.file, ' (', num2str(sweepNumber), ') - fft')); % naming figure file
% % % plot(f,P1) 
% % % title('Single-Sided Amplitude Spectrum of X(t)')
% % % xlabel('f (Hz)')
% % % ylabel('|P1(f)|')

% [yFiltered,d] = highpass(y,highpassThreshold,samplingFrequency);
% 
% [yupper,ylower] = envelope(yFiltered);
% figure
% plot(yFiltered)
% figure
% plot(yupper)

plot(y)
axis([50000 250000 -300 50]); 
yFiltered=highpass(y,highpassThreshold,samplingFrequency);
figure;
plot(yFiltered)
axis([50000 250000 -300 50]); 

yFilteredTwice=lowpass(yFiltered,lowpassThreshold,samplingFrequency);
figure;
plot(yFilteredTwice)
axis([50000 250000 -300 50]);
% for whatever reason, using the highpass filter first is better - the other
% way around creates more oscillations.

% lowpass(y,lowpassThreshold,samplingFrequency)
% bandpass(y,[highpassThreshold lowpassThreshold],samplingFrequency) %'ImpulseResponse','iir','Steepness', 0.9
% yFiltered = bandpass(y,[highpassThreshold lowpassThreshold],samplingFrequency);
% [yupper,ylower] = envelope(yFiltered);
% figure
% plot(yFiltered);
% figure
% plot(yupper);
% figure
% plot(abs(yFiltered))





end