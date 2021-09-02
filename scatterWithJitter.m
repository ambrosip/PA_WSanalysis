function scatterWithJitter(cohort, data)

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% organize data
scatterX = data(:,1);
scatterY = data(:,2);

%% generate fig 
figure('name', strcat("scatterWithJitter_", analysisDate, "_", cohort), 'Units', 'points');
hold on;
set(gcf,'OuterPosition',[0 100 220 222]);
set(gca,'FontName','Arial');
set(gca,'TickLength', [0.025, 0.025]);
set(gca,'LineWidth', 0.75);
axis([0.5 4.5 0 2500]); % for confocal data
% axis([0.5 4.5 0 5]);  % for oIPSC data

scatter(scatterX, scatterY, 25, [204/255, 121/255, 167/255], 'filled', 'MarkerFaceAlpha', 0.7, 'jitter', 'on', 'jitterAmount', 0.1);

ylabel('F/\mum^2','interpreter','Tex'); % for confocal data
% ylabel('oIPSC (nA)')  % for oIPSC data
xticks([0 1 2 3 4 5])
yticks([0:2500/3:2500]) % for confocal data
% yticks([0:5])   % for oIPSC data
xticklabels({'','a','a req','d','d req',''});
yticklabels({'0','','','','','5'});
set(findall(gcf,'-property','FontSize'),'FontSize',9)
hold off;

end