function scatterWithJitterDuo(cohort, data_not_connected, data_connected)

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% organize data
scatterX_notconnected = data_not_connected(:,1);
scatterY_notconnected = data_not_connected(:,2);
scatterX_connected = data_connected(:,1);
scatterY_connected = data_connected(:,2);

%% generate fig 
figure('name', strcat("scatterWithJitterDuo_", analysisDate, "_", cohort), 'Units', 'points');
hold on;
set(gcf,'OuterPosition',[0 100 220 222]);
set(gca,'FontName','Arial');
set(gca,'TickLength', [0.025, 0.025]);
set(gca,'LineWidth', 0.75);
% axis([0.5 4.5 0 2500]); % for confocal data
% axis([0.5 4.5 0 5]);    % for oIPSC data
% axis([0.5 4.5 -15 5]);    % for firing data SD median
axis([0.5 4.5 -25 5]);    % for firing data SD median

scatter(scatterX_connected, scatterY_connected, 25, [204/255, 121/255, 167/255], 'filled', 'MarkerFaceAlpha', 0.7, 'jitter', 'on', 'jitterAmount', 0.1);
scatter(scatterX_notconnected, scatterY_notconnected, 25, [0, 0, 0], 'filled', 'MarkerFaceAlpha', 0.7, 'jitter', 'on', 'jitterAmount', 0.1);

plot([0,5],[-2,-2],':', 'color', 'k');  % for SDs from firing data

xticks([0 1 2 3 4 5])
% xticklabels({'','dms','dls','','',''});
% xticklabels({'','a','a req','d','d req',''});
xticklabels({'','d req','d','dms','',''});
% xticklabels({'','a req','a','dls','',''});

% yticks([0:2500/3:2500]) % for confocal data
% yticks([0:5])   % for oIPSC data
yticks([-25:5:5])     % for Hz data

% yticklabels({'0','','','','','5'}); % for oIPSC data
% yticklabels({'-15','','','','5'}); % for Hz data
yticklabels({'-25','','','','','','5'}); % for Hz data

% ylabel('F/\mum^2','interpreter','Tex'); % for confocal data
% ylabel('oIPSC (nA)')    % for oIPSC data
% ylabel('SDs from avg baseline Hz')      % for Hz data
ylabel('SDs from pre-light Hz')      % for Hz data


set(findall(gcf,'-property','FontSize'),'FontSize',9)
hold off;

end