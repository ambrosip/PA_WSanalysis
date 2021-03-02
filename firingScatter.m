% INPUTS:
%   - cohort: a string to name the figure
%   - data: array with 4 columns properly ordered - I just saved the 4
%   columns from excel into a variable in the command window like this:
%       data = [*ctrl+V*] 

function firingScatter(cohort, data)

[preLightHzMean, preLightHzSD, duringLightHzMean, duringLightHzSD] = firingScatterPrep(data);

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% generate fig - vertical and horizontal error bars
figure('name', strcat("firingScatter1_", analysisDate, "_", cohort));
hold on;
errorbar(preLightHzMean, duringLightHzMean, duringLightHzSD, duringLightHzSD, preLightHzSD, preLightHzSD, 'o', 'color', 'k', 'MarkerFaceColor', 'k');
plot([0,15],[0,15],'--', 'color', 'k');
axis([0 15 0 15]);
set(gcf,'Position',[1000 50 400 400]);
hold off;

% generate fig - horizontal error bars only
figure('name', strcat("firingScatter2_", analysisDate, "_", cohort));
hold on;
errorbar(preLightHzMean, duringLightHzMean, preLightHzSD, 'horizontal', 'o', 'color', 'k', 'MarkerFaceColor', 'k');
plot([0,15],[0,15],'--', 'color', 'k');
axis([0 15 0 15]);
set(gcf,'Position',[500 50 400 400]);

% generate fig - horizontal error bars only
figure('name', strcat("firingScatter3_", analysisDate, "_", cohort), 'Units', 'points');
hold on;
set(gcf,'OuterPosition',[50 50 180 300]);
set(gca,'FontName','Arial');
set(gca,'TickLength', [0.025, 0.025]);
set(gca,'LineWidth', 0.75);
scatter(preLightHzMean, duringLightHzMean, 25, [0 0 0], 'filled', 'MarkerFaceAlpha', 0.7);
errorbar(preLightHzMean, duringLightHzMean, preLightHzSD, '.', 'horizontal', 'color', [0 0 0], 'CapSize', 0);
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
