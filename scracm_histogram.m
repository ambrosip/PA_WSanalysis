
close all
allDifferences = a;
ymaxhist = 100;
edges=-1:0.1:1;
figure('name', 'histogram_DMS_prob');
h=histogram(allDifferences,'BinEdges',edges,'Normalization','percentage','DisplayStyle', 'stairs', 'EdgeColor', 'k');
% h=histogram(allDifferences,10,'Normalization','percentage','DisplayStyle', 'stairs', 'EdgeColor', 'k');
% h=histogram('BinEdges', xMinInSec:xMaxInSec, 'BinCounts', firingHz, 'DisplayStyle', 'stairs', 'EdgeColor', 'k'); 
xlabel('P(oIPSC)');
ylabel('Percentage of ROIs');
% axis([-1 1 0 ymaxhist]);
axis([-inf inf 0 ymaxhist]);
% axis([0 1 0 ymaxhist]);
yticks([0 ymaxhist]);

% h.BinWidth=0.2
% h.BinLimits=[-1.1 1.1]


% h.BinWidth=0.1
h.BinLimits=[0 1]
                      