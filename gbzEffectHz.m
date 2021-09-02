function gbzEffectHz(cohort, notconnected_data, connected_data)

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');


%% generate fig with all data
figure('name', strcat("gbzEffectHz_", analysisDate, "_", cohort), 'Units', 'points');
hold on;
set(gcf,'OuterPosition',[0 100 110 222]);
set(gca,'FontName','Arial');
set(gca,'TickLength', [0.025, 0.025]);
set(gca,'LineWidth', 0.75);
axis([0.5 2.5 -0.1 1.2]);
plot([0.5 2.5],[1 1], ':', 'color', 'k');

% plot not connected cell data in black
for row = 1:size(notconnected_data,1)
    scatter([1:2], notconnected_data(row,:), 25, [0 0 0], 'filled', 'MarkerFaceAlpha', 0.7);
    plot(notconnected_data(row,:), 'color', [0 0 0 0.7]);
end

% plot connected cell data in pink
for row = 1:size(connected_data,1)
    scatter([1:2], connected_data(row,:), 25, [204/255, 121/255, 167/255], 'filled', 'MarkerFaceAlpha', 0.7);
    plot(connected_data(row,:), 'color', [204/255, 121/255, 167/255 0.7]);
end

ylabel('Hz ON / Hz OFF');
yticks([0:0.2:1.2])
yticklabels({'0','','','','','','1.2'});
xticklabels({'-GBZ','+GBZ'});
set(findall(gcf,'-property','FontSize'),'FontSize',9)
hold off;

end