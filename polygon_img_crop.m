%{ 
DOCUMENTATION
Created: 2022 09 13
Works? No
Author: PA

This function is used to crop an image from Ocular to fit the Polygon field.

INPUTS explained:

    - gridColumns: number of columns in the polygon ROI grid

    - gridRows: number of rows in the polygon ROI grid

    - cellImageFileName: file containing image you want to crop

    - cellImageDir: folder where cell image file is located

INPUTS defaults:      
    
OUTPUTS:
    Fig
    
ASSUMPTIONS: 

BEWARE:

TO DO:

%}

% function polygon_img_crop()
%%  USER INPUT ============================================================

gridColumns = 5;
gridRows = 5;
cellImageFileName = 'pre.tif';
cellImageDir = 'R:\Basic_Sciences\Phys\Lerner_Lab_tnl2633\Priscilla\Ephys\2022\20220913 polygon photobleaching\fifth test';


%% PREP - add path

addpath(cellImageDir);


%% PREP - get monitor info for plot display organization =====================================================

monitorPositions = get(0, 'MonitorPositions' );
maxHeight = monitorPositions(1,4) - 100;
maxWidth = monitorPositions(1,3) - 100;


%% PREP - get info from h5 file  ============================================================

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% % get info from h5 files in the directory
% % find all the H5 files in dir 
% files=dir(fullfile(dirName, '*.h5'));
% 
% % get the names of the files
% filesNames={files.name}';
% 
% % get mouse and data info from the first name
% fileName=filesNames{1};


%% PLOT 1 - cropped cell image ===============================================================================
% make sure you're getting the image taken with zoom = 1
% concatenate strings from user input to get full path to figure file
cellImageFileDir = strcat(cellImageDir,'\',cellImageFileName);

% import image
cellImage = imread(cellImageFileDir);

% crop image according to the polygon's mirror calibration
% I verified this in PPT lol
% the original image is 1376 x 1024 pixels
% ASSUMPTION ALERT: the calibration of the polygon will not change over time
% I need to crop the top 100 pixels and the bottom 51 pixels
% imcrop determines the rectangle to keep in the following form: [XMIN YMIN WIDTH HEIGHT]
% note that y increases from top to bottom, so ymin should match my
% required "top crop".
% I do not need to crop along the x axis, so xmin = 1 and width = 1376
% the height is 1024 - 100 - 51 = 873
croppedImage = imcrop(cellImage, [1,100,1376,873]);
figure('name', strcat(cellImageFileName, '_', analysisDate, ' - cropped image'));
hold on;
imshow(croppedImage, 'Border', 'tight');

% calculate parameters for scale bar
% ASSUMPTION ALERT: pixelsPerMicron corresponds to my usual 40x objective at 1x zoom
xmaxImg = size(croppedImage,2);    % in pixels, should be 1376
ymaxImg = size(croppedImage,1);    % in pixels, should be 874
scaleDistanceFromEdge = 50;     % in pixels
scaleBarSize = 50;              % in um
pixelsPerMicron = 873 / 222.2;  % 222.2 um in 873 pixels
scaleBarSizeInPixels = scaleBarSize * pixelsPerMicron;

% add scale bar to figure
% line([x x], [y y])
line([scaleDistanceFromEdge scaleDistanceFromEdge+scaleBarSizeInPixels],...
    [ymaxImg-scaleDistanceFromEdge ymaxImg-scaleDistanceFromEdge],...
    'LineWidth', 2, 'Color', 'k');
% text(x, y, string)
text(scaleDistanceFromEdge,...
    ymaxImg-2*scaleDistanceFromEdge,...
    strcat(num2str(scaleBarSize), " μm"),...
    'FontSize', 10);
hold off;

% get inner figure size and store half of those values
pos = get(gcf, 'InnerPosition'); %// gives x left, y bottom, width, height
innerWidth = pos(3)/2;
innerHeight = pos(4)/2;

% get outer figure size and store half of those values
pos = get(gcf, 'OuterPosition'); %// gives x left, y bottom, width, height
outerWidth = pos(3)/2;
outerHeight = pos(4)/2;

% set figure size to the stored values
set(gcf,'InnerPosition',[0 maxHeight-innerHeight innerWidth innerHeight]);


%% PLOT 2 - cropped cell image with grid ===============================================================================
% make sure you're getting the image taken with zoom = 1
% concatenate strings from user input to get full path to figure file
cellImageFileDir = strcat(cellImageDir,'\',cellImageFileName);

% import image
cellImage = imread(cellImageFileDir);

% crop image according to the polygon's mirror calibration
croppedImage = imcrop(cellImage, [1,100,1376,873]);
figure('name', strcat(cellImageFileName, '_', analysisDate, ' - cropped image grid'));
hold on;
imshow(croppedImage, 'Border', 'tight');

% add polygon grid to figure
for row = 1:gridRows+1
    ypos = (row-1)*ymaxImg/gridRows;
    yline(ypos, 'Color', 'k', 'LineWidth', 0.5);
end
for col = 1:gridColumns+1
    xpos = (col-1)*xmaxImg/gridColumns;
    xline(xpos, 'Color', 'k', 'LineWidth', 0.5);
end

% add scale bar to figure
% line([x x], [y y])
line([scaleDistanceFromEdge scaleDistanceFromEdge+scaleBarSizeInPixels],...
    [ymaxImg-scaleDistanceFromEdge ymaxImg-scaleDistanceFromEdge],...
    'LineWidth', 2, 'Color', 'k');
% text(x, y, string)
text(scaleDistanceFromEdge,...
    ymaxImg-2*scaleDistanceFromEdge,...
    strcat(num2str(scaleBarSize), " μm"),...
    'FontSize', 10, 'Color', 'k');
hold off;

% re-size 
set(gcf,'InnerPosition',[0 maxHeight-innerHeight innerWidth innerHeight]);

