% INPUTS:
%   - cohort: a string to name the figure
%   - data: array with 4 columns properly ordered - I just saved the 4
%   columns from excel into a variable in the command window like this:
%       data = [*ctrl+V*] 

function firingScatterDuo(cohort, notconnected_data, connected_data)

[notconnected_preLightHzMean, notconnected_preLightHzSD, notconnected_duringLightHzMean, notconnected_duringLightHzSD] = firingScatterPrep(notconnected_data);
[connected_preLightHzMean, connected_preLightHzSD, connected_duringLightHzMean, connected_duringLightHzSD] = firingScatterPrep(connected_data);

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% % generate fig - vertical and horizontal error bars
% figure('name', strcat("firingScatterDuo_1_", analysisDate, "_", cohort), 'Units', 'points');
% hold on;
% set(gcf,'OuterPosition',[1000 50 180 300]);
% set(gca,'FontName','Arial');
% errorbar(notconnected_preLightHzMean, notconnected_duringLightHzMean, notconnected_duringLightHzSD, notconnected_duringLightHzSD, notconnected_preLightHzSD, notconnected_preLightHzSD, 'o', 'color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 4);
% errorbar(connected_preLightHzMean, connected_duringLightHzMean, connected_duringLightHzSD, connected_duringLightHzSD, connected_preLightHzSD, connected_preLightHzSD, 'o', 'color', [204/255, 121/255, 167/255], 'MarkerFaceColor', [204/255, 121/255, 167/255], 'MarkerSize', 4);
% plot([0,15],[0,15],'--', 'color', 'k');
% axis([0 15 0 15]);
% xlabel('Hz before o-stim');
% ylabel('Hz during o-stim');
% xticklabels({'0','','','15'});
% yticklabels({'0','','','15'});
% set(findall(gcf,'-property','FontSize'),'FontSize',9)
% axis square
% box on;
% hold off;

% % generate fig - horizontal error bars only
% figure('name', strcat("firingScatterDuo_2_", analysisDate, "_", cohort), 'Units', 'points');
% hold on;
% set(gcf,'OuterPosition',[500 50 180 300]);
% set(gca,'FontName','Arial');
% errorbar(notconnected_preLightHzMean, notconnected_duringLightHzMean, notconnected_preLightHzSD, 'horizontal', 'o', 'color', [0 0 0], 'MarkerFaceColor', [0 0 0], 'MarkerSize', 4);
% errorbar(connected_preLightHzMean, connected_duringLightHzMean, connected_preLightHzSD, 'horizontal', 'o', 'color', [204/255, 121/255, 167/255], 'MarkerFaceColor', [204/255, 121/255, 167/255], 'MarkerSize', 4);
% plot([0,15],[0,15],'--', 'color', 'k');
% axis([0 15 0 15]);
% xlabel('Hz before o-stim');
% ylabel('Hz during o-stim');
% xticklabels({'0','','','15'});
% yticklabels({'0','','','15'});
% set(findall(gcf,'-property','FontSize'),'FontSize',9)
% axis square
% box on;
% hold off;

% % generate fig - horizontal error bars only
% figure('name', strcat("firingScatterDuo_3_", analysisDate, "_", cohort), 'Units', 'points');
% hold on;
% set(gcf,'OuterPosition',[50 50 180 300]);
% set(gca,'FontName','Arial');
% set(gca,'TickLength', [0.025, 0.025]);
% set(gca,'LineWidth', 0.75);
% scatter(notconnected_preLightHzMean, notconnected_duringLightHzMean, 25, [0 0 0], 'filled', 'MarkerFaceAlpha', 0.7);
% scatter(connected_preLightHzMean, connected_duringLightHzMean, 25, [204/255, 121/255, 167/255], 'filled', 'MarkerFaceAlpha', 0.7);
% errorbar(notconnected_preLightHzMean, notconnected_duringLightHzMean, notconnected_preLightHzSD, '.', 'horizontal', 'color', [0 0 0], 'CapSize', 0);
% errorbar(connected_preLightHzMean, connected_duringLightHzMean, connected_preLightHzSD, '.', 'horizontal', 'color', [204/255, 121/255, 167/255], 'CapSize', 0);
% plot([0,15],[0,15],':', 'color', 'k');
% plot([-1,0],[-1,-1],'color',[1 1 1],'LineWidth',2);
% plot([-1,-1],[-1,0],'color',[1 1 1],'LineWidth',2);
% axis([-1 15 -1 15]);
% xlabel('Hz before o-stim');
% ylabel('Hz during o-stim');
% xticklabels({'0','','','15'});
% yticklabels({'0','','','15'});
% set(findall(gcf,'-property','FontSize'),'FontSize',9)
% axis square
% % box on;
% hold off;

% generate fig - horizontal and vertical error bars
figure('name', strcat("firingScatterDuo_4_", analysisDate, "_", cohort), 'Units', 'points');
hold on;
set(gcf,'OuterPosition',[50 50 180 300]);
set(gca,'FontName','Arial');
set(gca,'TickLength', [0.025, 0.025]);
set(gca,'LineWidth', 0.75);
% dots
scatter(notconnected_preLightHzMean, notconnected_duringLightHzMean, 25, [0 0 0], 'filled', 'MarkerFaceAlpha', 0.7);
scatter(connected_preLightHzMean, connected_duringLightHzMean, 25, [204/255, 121/255, 167/255], 'filled', 'MarkerFaceAlpha', 0.7);
% horizontal error
errorbar(notconnected_preLightHzMean, notconnected_duringLightHzMean, notconnected_preLightHzSD, '.', 'horizontal', 'color', [0 0 0], 'CapSize', 0);
errorbar(connected_preLightHzMean, connected_duringLightHzMean, connected_preLightHzSD, '.', 'horizontal', 'color', [204/255, 121/255, 167/255], 'CapSize', 0);
% vertical error
errorbar(notconnected_preLightHzMean, notconnected_duringLightHzMean, notconnected_duringLightHzSD, '.', 'vertical', 'color', [0 0 0], 'CapSize', 0);
errorbar(connected_preLightHzMean, connected_duringLightHzMean, connected_duringLightHzSD, '.', 'vertical', 'color', [204/255, 121/255, 167/255], 'CapSize', 0);
plot([0,15],[0,15],':', 'color', 'k');
plot([-1,0],[-1,-1],'color',[1 1 1],'LineWidth',2);
plot([-1,-1],[-1,0],'color',[1 1 1],'LineWidth',2);
axis([-1 15 -1 15]);
xlabel('Hz before o-stim');
ylabel('Hz during o-stim');
xticklabels({'0','','','15'});
yticklabels({'0','','','15'});
set(findall(gcf,'-property','FontSize'),'FontSize',9)
axis square
% box on;
hold off;

end
