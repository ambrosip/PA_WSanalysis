%{ 
DOCUMENTATION
Created: 2024 10 28
Works? Yes
Author: PA

This function is used to crop an image from Ocular to fit the Polygon field.

FUNCTION INPUTS explained:

    - cellImageFile: file containing image you want to crop. Don't add
    ".tif" extension - the code will add it for you. Use quotes!

    - cellImageDir: folder where cell image file is located. Use quotes!

    - leftCrop: the polygon mirror area is smaller than the field of view
    of the rig camera. To properly overlay the grid onto the images taken
    from the rig, I need to crop the image from ocular. This crop requires
    four parameters: leftCrop, rightCrop, topCrop and bottomCrop. The first
    corresponds to the number of pixels in the x-axis betwween the left-most 
    edge of the photo and the left-most edge of the mirror grid.

    - rightCrop: number of pixels between the right-most edge of the
    polygon's mirror grid and the right edge of the ocular picture.

    - topCrop: number of pixels between the top-most edge of the
    polygon's mirror grid and the top edge of the ocular picture.

    - bottomCrop: number of pixels between the bottom edge of the
    polygon's mirror grid and the bottom edge of the ocular picture.

    - gridColumns: number of columns in the polygon ROI grid

    - gridRows: number of rows in the polygon ROI grid


USER INPUTS explained:

INPUTS defaults:      
    
OUTPUTS:
    Fig
    
ASSUMPTIONS: 

BEWARE:

TO DO:
%}



function polygon_img_crop(cellImageFile, cellImageDir, leftCrop, rightCrop, topCrop, bottomCrop, gridColumns, gridRows, prefix, renameFiles, cellNamePrefix)
%%  USER INPUT ============================================================

gridFillHorizontal = 1;     % 0.067 for 15 column
gridFillVertical = 1;       % 0.067 for 15 line

% saveDir = '/Users/priscilla/OHSU Dropbox/Priscilla Ambrosi/Dropbox - Lerner Lab/Ambrosi et al_sCRACM 2024/Data Analysis/2024-10-28';


%% PREP - add path and prep cellImageFileName ===========================================================================================

addpath(cellImageDir);

% add tif to end of cellImageFile
cellImageFileName = strcat(cellImageFile, '.tif');


%% PREP - get monitor info for plot display organization =====================================================

monitorPositions = get(0, 'MonitorPositions' );
maxHeight = monitorPositions(1,4) - 100;
maxWidth = monitorPositions(1,3) - 100;


%% PREP - get info from h5 file  ============================================================

% get today's date for naming output files
analysisDate =  datestr(datetime('today'),'yyyy-mm-dd');

% get info from h5 files in the directory
% find all the H5 files in dir 
files=dir(fullfile(cellImageDir, '*.h5'));

% get the names of the files
filesNames={files.name}';

% get mouse and data info from the first name
fileName=filesNames{1};
fileName=fileName(1,1:16);


%% PLOT 1 - cropped image with white grid ===============================================================================
% make sure you're getting the image taken with zoom = 1
% concatenate strings from user input to get full path to figure file
% add if statement to add OS-appropriate slash
if ismac
    cellImageFileDir = strcat(cellImageDir,'/',cellImageFileName);
elseif ispc
    cellImageFileDir = strcat(cellImageDir,'\',cellImageFileName);
end

% import image
cellImage = imread(cellImageFileDir);

% crop image according to the polygon's mirror calibration
% Old code (prior to 2022 12 21):
% % I verified this in PPT lol
% % the original image is 1376 x 1024 pixels
% % ASSUMPTION ALERT: the calibration of the polygon will not change over time
% % I need to crop the top 100 pixels and the bottom 51 pixels
% % imcrop determines the rectangle to keep in the following form: [XMIN YMIN WIDTH HEIGHT]
% % note that y increases from top to bottom, so ymin should match my
% % required "top crop".
% % I do not need to crop along the x axis, so xmin = 1 and width = 1376
% % the height is 1024 - 100 - 51 = 873
% croppedImage = imcrop(cellImage, [1,100,1376,873]);

% custom-crop based on user input (in case the polygon's mirror calibration changes,
% like it did on 2022 12 21)
cropXmin = 1 + leftCrop;
cropWidth = 1376 - rightCrop - leftCrop;
cropYmin = 1 + topCrop;
cropHeight = 1024 - topCrop - bottomCrop;
croppedImage = imcrop(cellImage, [cropXmin, cropYmin, cropWidth, cropHeight]);

% name and create figure
if renameFiles == 1
    fig1 = figure('name', strcat(cellNamePrefix, prefix));
else
    fig1 = figure('name', strcat(fileName, prefix, cellImageFileName, '_cropped_on_', analysisDate));
end
ax1 = axes(fig1);
hold on;
imshow(croppedImage, 'Border', 'tight');

% calculate parameters for scale bar
% ASSUMPTION ALERT: pixelsPerMicron corresponds to my usual 40x objective at 1x zoom
xmaxImg = size(croppedImage,2);    % in pixels, should be 1376 (prior to 2022 12 21)
ymaxImg = size(croppedImage,1);    % in pixels, should be 874 (prior to 2022 12 21)
scaleDistanceFromEdge = 50;     % in pixels
scaleBarSize = 50;              % in um
pixelsPerMicron = 873 / 222.2;  % 222.2 um in 873 pixels
scaleBarSizeInPixels = scaleBarSize * pixelsPerMicron;

% add polygon grid to figure
% also keep track of points to make rectangles if fill is NOT 1 (see below)
allYpos = [];
for row = 1:gridRows+1
    ypos = (row-1)*ymaxImg/gridRows;
    allYpos = [allYpos; ypos];
    if gridFillHorizontal == 1 && gridFillVertical == 1
        yline(ypos, 'Color', 'w', 'LineWidth', 0.5);
    end
end

allXpos = [];
for col = 1:gridColumns+1
    xpos = (col-1)*xmaxImg/gridColumns;
    allXpos = [allXpos; xpos];
    if gridFillHorizontal == 1 && gridFillVertical == 1
        xline(xpos, 'Color', 'w', 'LineWidth', 0.5);
    end
end

% add rectangles to figure if the grid fill is NOT 1 (aka the o-stim ROI is
% smaller than the full ROI)
% first figure out the dimensions of the rectangle
ROIwidth = gridFillHorizontal * xmaxImg/gridColumns;
ROIheight = gridFillVertical * ymaxImg/gridRows;
extraWidth = (1 - gridFillHorizontal) * xmaxImg/gridColumns;
extraHeight = (1 - gridFillVertical) * ymaxImg/gridRows;
if gridFillHorizontal ~= 1 || gridFillVertical ~= 1
    for row = 1:length(allYpos)-1
        for col = 1:length(allXpos)-1
            ROIx = allXpos(col) + extraWidth/2;
            ROIy = allYpos(row) + extraHeight/2;
            rectangle('Position', [ROIx, ROIy, ROIwidth, ROIheight]);
        end
    end
end

% add scale bar to figure
% line([x x], [y y])
line([scaleDistanceFromEdge scaleDistanceFromEdge+scaleBarSizeInPixels],...
    [ymaxImg-scaleDistanceFromEdge ymaxImg-scaleDistanceFromEdge],...
    'LineWidth', 2, 'Color', 'w');
% text(x, y, string)
text(scaleDistanceFromEdge,...
    ymaxImg-2*scaleDistanceFromEdge,...
    strcat(num2str(scaleBarSize), " μm"),...
    'FontSize', 10,'Color','w');
hold off;

% % get inner figure size and store half of those values
% pos = get(gcf, 'InnerPosition'); %// gives x left, y bottom, width, height
% innerWidth = pos(3)/2;
% innerHeight = pos(4)/2;
% 
% % get outer figure size and store half of those values
% pos = get(gcf, 'OuterPosition'); %// gives x left, y bottom, width, height
% outerWidth = pos(3)/2;
% outerHeight = pos(4)/2;
% 
% % set figure size to the stored values
% set(gcf,'InnerPosition',[0 maxHeight-innerHeight innerWidth innerHeight]);


%% PLOT 2 - cropped inverted image with grid ===============================================================================
% normalize image to max intensity value
% this is important for dealing with summed z-stacks (instead of max
% intensity z-stacks)
if max(max(cellImage))/256 > 1
    adjustmentFactor = (max(max(cellImage))/256) * 256;
    cellImage = cellImage/adjustmentFactor;
end

% invert image so that black = white
% I wanted to do this anyway, but MATLAB also forced my hand. When I save
% images with my saveAllFigs function, MATLAB turns all white annotations
% (text and lines) from white to black, rendering my scale bar useless in
% tif files. The lines are still there in svg files!
invertedImage = imcomplement(cellImage);

% crop image according to the polygon's mirror calibration
croppedImage = imcrop(invertedImage, [cropXmin, cropYmin, cropWidth, cropHeight]);

% create fig
if renameFiles == 1
    fig2 = figure('name', strcat(cellNamePrefix, prefix, ' - inverted'));
else
    fig2 = figure('name', strcat(fileName, prefix, cellImageFileName, '_cropped_on_', analysisDate, ' - inverted'));
end
ax2 = axes(fig2);
hold on;
imshow(croppedImage, 'Border', 'tight');

% add polygon grid to figure
% also keep track of points to make rectangles if fill is NOT 1 (see below)
allYpos = [];
for row = 1:gridRows+1
    ypos = (row-1)*ymaxImg/gridRows;
    allYpos = [allYpos; ypos];
    if gridFillHorizontal == 1 && gridFillVertical == 1
        yline(ax2, ypos, 'Color', 'k', 'LineWidth', 0.5);
    end
end

allXpos = [];
for col = 1:gridColumns+1
    xpos = (col-1)*xmaxImg/gridColumns;
    allXpos = [allXpos; xpos];
    if gridFillHorizontal == 1 && gridFillVertical == 1
        xline(ax2, xpos, 'Color', 'k', 'LineWidth', 0.5);
    end
end

% add rectangles to figure if the grid fill is NOT 1 (aka the o-stim ROI is
% smaller than the full ROI)
if gridFillHorizontal ~= 1 || gridFillVertical ~= 1
    for row = 1:length(allYpos)-1
        for col = 1:length(allXpos)-1
            ROIx = allXpos(col) + extraWidth/2;
            ROIy = allYpos(row) + extraHeight/2;
            rectangle(ax2, 'Position', [ROIx, ROIy, ROIwidth, ROIheight]);
        end
    end
end

% add scale bar to figure
% line([x x], [y y])
line(ax2, [scaleDistanceFromEdge scaleDistanceFromEdge+scaleBarSizeInPixels],...
    [ymaxImg-scaleDistanceFromEdge ymaxImg-scaleDistanceFromEdge],...
    'LineWidth', 2, 'Color', 'k');
% text(x, y, string)
text(ax2, scaleDistanceFromEdge,...
    ymaxImg-2*scaleDistanceFromEdge,...
    strcat(num2str(scaleBarSize), " μm"),...
    'FontSize', 10,'Color','k');
hold off;


%% SAVE figs ==============================================================

% saveAllFigs(saveDir);