function [allRs, allSweeps] = rs(objectArray,varargin)

    % optional arguments
    % set defaults for optional inputs 
%     optargs = {'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries\From MATLAB'};
    optargs = {'D:\Temp\From MATLAB'};
    % overwrite defaults with values specified in varargin
    numvarargs = length(varargin);
    optargs(1:numvarargs) = varargin;
    % place optional args in memorable variable names
    [savefileto] = optargs{:};

% ex: rs([m012.s0001; m012.s0004])
% ex: rs([m012.s0001; m012.s0004],'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Data summaries')

% uncomment below if you want to look at the test pulse raw data
% arrayfun(@(obj) plotxyallch(obj,-inf, inf, -1000, 800, false), objectArray);

[allRs, allSweeps] = arrayfun(@(obj) calculateRs(obj), objectArray);

figure('name', strcat(objectArray(1).file(1:15)," - Rs for ",num2str(allSweeps(1)),'-',num2str(allSweeps(end)))); % naming figure file
plot(allSweeps, allRs,'-o');

% plot lines marking 30% increase and 30% decrese in Rs compared to first
% test pulse
line([allSweeps(1)-5 allSweeps(end)+5],[allRs(1)*0.7, allRs(1)*0.7],'Color','black','LineStyle','--')
line([allSweeps(1)-5 allSweeps(end)+5],[allRs(1)*1.3, allRs(1)*1.3],'Color','black','LineStyle','--')
axis([allSweeps(1)-5 allSweeps(end)+5 0 60])
ylabel('Rs (M\Omega)');
xlabel('Sweeps');
title([strcat(objectArray(1).file(1:15)," - Rs: ",num2str(allSweeps(1)),'-',num2str(allSweeps(end)))],'Interpreter','none');

% % save csv file with sweep number and Rs 
% % this works but does not let me label the data
% if nargin>1
%     savefileto = varargin;
%     filename = strcat(objectArray(1).file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - Rs");
%     fulldirectory = strcat(savefileto,'\',filename,'.csv');
%     csvwrite(fulldirectory,[allSweeps;allRs]);
% else
%     disp('gimme directory if you want this saved')
% end

    % save csv file with sweep number and Rs 
        filename = strcat(objectArray(1).file(1:15),'_',num2str(allSweeps(1)),'-',num2str(allSweeps(end))," - Rs");
        fulldirectory = strcat(savefileto,'\',filename,'.csv');
        dataInCellFormat = {};
        dataInCellFormat = num2cell([transpose(allSweeps), transpose(allRs)]);
        labeledData = cell2table(dataInCellFormat,'VariableNames',{'sweep','rs'});
        writetable(labeledData,fulldirectory);     
        disp('I saved it, ur welcome love')
        disp('Change directory if you want this saved elsewhere!')


end