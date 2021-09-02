function ipscScatter(cohort, data)

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

ipsc_latency = data(:,1);
ipsc_amplitude = data(:,2);

% generate fig - vertical and horizontal error bars
figure('name', strcat("ipscScatter_", analysisDate, "_", cohort), 'Units', 'points');
hold on;
set(gcf,'OuterPosition',[50 50 110 222]);
set(gca,'FontName','Arial');
set(gca,'TickLength', [0.025, 0.025]);
set(gca,'LineWidth', 0.75);
scatter(ipsc_latency, ipsc_amplitude, 25, [204/255, 121/255, 167/255], 'filled', 'MarkerFaceAlpha', 0.7);
plot([1,1],[0,5],':', 'color', 'k');
axis([-0.25 5 -0.5 5]);
% xlabel('Latency (ms)');
ylabel('oIPSC (nA)');
xticks([0 1 2 3 4 5])
yticks([0 1 2 3 4 5])
xticklabels({'0','','','','','5'});
yticklabels({'0','','','','','5'});
set(findall(gcf,'-property','FontSize'),'FontSize',9)
hold off;

end