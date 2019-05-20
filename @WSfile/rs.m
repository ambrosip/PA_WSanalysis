function rs(objectArray)

% ex: rs([m012.s0001; m012.s0004])

% uncomment if you want to look at the test pulse raw data
arrayfun(@(obj) plotxyallch(obj,-inf, inf, -1000, 800, false), objectArray);

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

end