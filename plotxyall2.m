function plotxyall2(mouseNumber, firstSweep, varargin)
   
    strcat('m', num2str(mouseNumber, '%03.f')) = loadedMouse
    loadedMouse = struct
    loadedMouse = strcat('m', num2str(mouseNumber, '%03.f'))
    loadedMouse
    sNumber = strcat('s', num2str(firstSweep, '%04.f'))
    loadedS = struct
    loadedS = getfield(loadedMouse, sNumber)
    fileName = getfield(loadedS, file)
    
    
%     
%     
%     strcat('s', num2str(firstSweep, '%04.f')) = struct
%     loadedS = struct;
%     loadedS = strcat('s', num2str(firstSweep, '%04.f'))
%     
%     fileName = getfield(loadedS,
%     
%     
%     loadedSweep = struct;
%     loadedSweep = getfield(sNumber,
%     loadedSweep = loadedMouse.(sNumber)
%     
%      
%      loadedMouse = struct;
%      loadedS = struct;
%      loadedSweep = struct;
% 
%      mouseName = strcat('m', num2str(mouseNumber, '%03.f'));
%      sNumber = strcat('s', num2str(firstSweep, '%04.f'));
%      sweepFirstName = strcat('sweep_', num2str(firstSweep, '%04.f'));
%      
%      varName = matlab.lang.makeValidName(mouseName)
%      
%      mouseName
%      sNumber
%      sweepFirstName
%      
% %      loadedMouse = getfield(mouseName, sNumber);
% %      loadedS = getfield(loadedMouse, sNumber);
% %      loadedMouse.(loadedS) = mouseName('s0001');
% %      loadedS = sNumber;
%      loadedMouse.(loadedS) = mouseName.(sNumber);
%      loadedMouse
%      fileName = getfield(loadedS, file);
%      
%      if length(fileName) == 28
%          firstSweepNumber = str2num(fileName(end-11:end-8));
%          lastSweepNumber = str2num(fileName(end-6:end-3));
%          sweepLastName = strcat('sweep_', num2str(lastSweepNumber));
%          
%      elseif length(fileName) < 28
%          firstSweepNumber = str2num(fileName(end-6:end-3));
%          lastSweepNumber = str2num(mfileName(end-6:end-3));
%          sweepLastName = strcat('sweep_', num2str(lastSweepNumber));
%   
%      end
%       
% %      samplingFrequency = mouseName.sNumber.header.Acquisition.SampleRate;
%      samplingFrequency = loadedMouse.loadedS.header.Acquisition.SampleRate;
%      
%      currentSweep = firstSweepNumber;
%      while currentSweep <= lastSweepNumber
%          currentSweepName = strcat('sweep_', num2str(currentSweep));
%          
%          loadedSweep = getfield(loadedS, currentSweepName);
%          y = loadedSweep.analogScans(:,1);
%          sweepDuration = size(y,1)/samplingFrequency;
%          x = linspace(0,sweepDuration,size(y,1))';
%          figure;
%          plot(x,y);
% 
%          numvarargs = length(varargin);
%          optargs = {-inf inf -inf inf};
%          optargs(1:numvarargs) = varargin;
%          [xmin, xmax, ymin, ymax] = optargs{:};
%          axis([xmin xmax ymin ymax])
%          title([fileName ' (' num2str(currentSweep) ')'],'Interpreter','none');
%             
%          currentSweep = currentSweep + 1;
     end
     
      