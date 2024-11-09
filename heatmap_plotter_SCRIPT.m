%% Heatmap plotter

%% COLOR MAPS

% make a color map from white to pink (reddish purple used in my paper)
% for coloring success rate heatmap
customColorMapPink = [linspace(1,204/255,256)' linspace(1,121/255,256)'  linspace(1,167/255,256)'];

% make color map from white to vermillion
% for coloring amplitude heatmap
customColorMapVermillion = [linspace(1,213/255,256)' linspace(1,94/255,256)'  linspace(1,0/255,256)'];

% make color map from white to blue
% for coloring onset latency heatmap
customColorMapBlue = [linspace(1,0/255,256)' linspace(1,114/255,256)'  linspace(1,178/255,256)'];

% make color map from white to sky blue
% for coloring onset latency jitter heatmap
customColorMapSky = [linspace(1,86/255,256)' linspace(1,180/255,256)'  linspace(1,233/255,256)'];


%% USER INPUT

fileName = "DLS charge 11";
gridColumns = 5;

data = [-1.395583831
1.462061013
20.16957327
269.1716853
571.8620352
3.901190044
2.672372396
247.581724
730.3691447
642.2853708
22.83733762
224.9259151
658.2883853
530.884407
46.58557992
65.03587502
436.0195361
118.1861431
8.358286257
7.263175734
54.35374747
62.72955424
6.028239292
12.12763765
4.896397561]';
% data = [0.32	0.62	0.64	0.76	0.46	0.48	0.58	0.62	0.56	0.28	0.36	0.6	0.8	0.42	0.26	0.54	0.66	0.78	0.34	0.28	0.48	0.5	0.56	0.46	0.26];
colorMap = customColorMapVermillion;

% set hetmap edges
heatmapMin = 0;
heatmapMax = 1800;


%% PREP

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% get info about heatmap size
leftCrop = 0;     % old: 0    old (2022-12-25): 0     new (2023-08-20): 0    % in pixels        % ALERT! this is new (2023-08-20)   
rightCrop = 14;    % old: 0    old (2022-12-25): 14    new (2023-08-20): 8    % in pixels        % ALERT! this is new (2023-08-20)
topCrop = 106;     % old: 100  old (2022-12-25): 106   new (2023-08-20): 69   % in pixels        % ALERT! this is new (2023-08-20)
bottomCrop = 58;  % old: 51   old (2022-12-25): 58    new (2023-08-20): 88   % in pixels        % ALERT! this is new (2023-08-20)

% get info about heatmap size
cropXmin = 1 + leftCrop;
cropWidth = 1376 - rightCrop - leftCrop;
cropYmin = 1 + topCrop;
cropHeight = 1024 - topCrop - bottomCrop;
% croppedImage = imcrop(cellImage, [cropXmin, cropYmin, cropWidth, cropHeight]);

% get monitor info
monitorPositions = get(0, 'MonitorPositions' );
maxHeight = monitorPositions(1,4) - 100;
maxWidth = monitorPositions(1,3) - 100;

% % get inner figure size and store half of those values
% pos = get(gcf, 'InnerPosition'); %// gives x left, y bottom, width, height
% innerWidth = pos(3)/2;
% innerHeight = pos(4)/2;


%% PLOT 9.2 - heatmap of average charge (AUC) normalized to largest charge

% organize data for heatmap
dataForHeatmap = reshape(data,gridColumns,[]).';

% resize heatmap
resizedHeatmap = imresize(dataForHeatmap, [cropHeight cropWidth], 'nearest');

% make heatmap without the heatmap function
% figure('name', strcat(fileName, " ", analysisDate, ' - avg charge heatmap')); % naming figure file
figure('name', strcat(fileName, " ", analysisDate, ' - charge')); % naming figure file
imshow(resizedHeatmap,'Colormap',colorMap,'DisplayRange', [heatmapMin,heatmapMax], 'Border', 'tight');
% title('ALl cells normalized avg charge')
title('charge')
% set(gcf,'InnerPosition',[innerWidth maxHeight-innerHeight innerWidth innerHeight]);
c = colorbar;
% c.Label.String = 'Normalized AVG charge';
c.Label.String = 'pA';