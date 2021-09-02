function gbzEffectIPSC(cohort, data)

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');


%% generate fig with all data
figure('name', strcat("gbzEffectIPSC_", analysisDate, "_", cohort), 'Units', 'points');
hold on;
set(gcf,'OuterPosition',[0 100 110 222]);
set(gca,'FontName','Arial');
set(gca,'TickLength', [0.025, 0.025]);
set(gca,'LineWidth', 0.75);
axis([0.5 2.5 -0.5 5]);

for row = 1:size(data,1)
    scatter([1:2], data(row,:), 25, [204/255, 121/255, 167/255], 'filled', 'MarkerFaceAlpha', 0.7);
    plot(data(row,:), 'color', [204/255, 121/255, 167/255 0.7]);
end

ylabel('oIPSC (nA)');
yticks([0 1 2 3 4 5])
yticklabels({'0','','','','','5'});
xticklabels({'-GBZ','+GBZ'});
set(findall(gcf,'-property','FontSize'),'FontSize',9)
hold off;

end