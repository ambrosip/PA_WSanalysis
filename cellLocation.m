% Use this code to illustrate the position of recorded cells in the slice.
% This code depends on a consistent naming strategy for all files.
% Slice images need to be saved as "s1_dic.tif"
% Cell location images need to be saved as "s1c1_5x_dic.tif" or "s1c1_4x_dic.tif"
% Wavesurfer files (.h5) need to be saved as "m027_2020-09-23_0001-0003.h5"
% FYI the code will fail if you have more than 9 slices imaged - but if you
% ever get to that point, you need to chill on your patching.
% Visual inspection is required to check the quality of the code output.

% INPUTS: 
%   - directory to be analyzed (dirName) when calling function
%       Ex: cellLocation('D:\Priscilla - BACKUP 20200319\Ephys\2020\20200731')
%   - coordinates of pipette location 
%   - registration parameters

% OUTPUTS:
%   - xlsx file with registered pipette locations
%   - figures
%       - raw slice images
%       - raw cell images + pipette location based on user input (red
%         circle)
%       - cell images registered to slice image (for quality control)
%       - slice image + cell locations (numbered white circles)

% TO DO: 
%   - test if any images need to be reflected prior to registration? So
%     far I am manually checking and reflecting images that need to be
%     fixed.


function cellLocation(dirName)

addpath(dirName);

%% USER INPUT

% (x,y) coordinates in pixel # of approximate pipette tip location
% use Fiji/ImageJ to estimate these coordinates.
pipx = 670; % 650 out of 1376 pixels
pipy = 526; % 515 out of 1024 pixels

% adjust registration parameters - monomodal 
% play with the parameters to get a good balance between precision and
% computation time. The columns on the right are parameters that worked
% well for some subsets of images.
[optimizer, metric] = imregconfig('monomodal');
optimizer.GradientMagnitudeTolerance = 1e-20;   % 1e-4      1e-10    1e-5
optimizer.MinimumStepLength = 0.1 ;             % 1e-5      1e-4     0.1
optimizer.MaximumStepLength = 1;             % 0.0625    0.06       1
optimizer.MaximumIterations = 10000;              % 100       500   10000
optimizer.RelaxationFactor = 0.5;               % 0.5       0.7     0.7

% adjust registration parameters - multimodal
% I started using multimodal but monomodal yields better results - so I
% commented out the multimodal code.
% [optimizer, metric] = imregconfig('multimodal');
% optimizer.InitialRadius = 0.002;    % 0.009     0.005   0.002
% optimizer.Epsilon = 1.5e-8;         % 1.5e-4    1.5e-6  1.5e-6 
% optimizer.GrowthFactor = 1.01;      % 1.05      1.01    1.01
% optimizer.MaximumIterations = 500;  % 100       300     300


%% MAIN CODE

% Matrices that will be filled
sliceImgsNames = [];
sliceImgs = [];
cellImgsNames = [];
cellImgs = [];
data = [];

% Get mouse number and date of experiment based on *.h5 files in folder
h5file = dir(fullfile(dirName, '*.h5')).name;
mouseNumber = h5file(1:4);
expDate = h5file(6:15);

% Convert strings to numbers so we can save data as xlsx later
expDateAsNum = str2num(strcat(h5file(6:9),h5file(11:12),h5file(14:15)));
mouseNumberAsNum = str2num(h5file(2:4));

% Get a list of all files in dir with the desired file name pattern.
% First, get slice images
% "?" is a wildcard/placeholder for a single character
% sliceFilePattern = fullfile(dirName, 's?_DLS_dic.tif');     % MODIFIED for m006 and m005
% sliceFilePattern = fullfile(dirName, 's?_DMS_dic.tif');     % MODIFIED for m006 and m005
sliceFilePattern = fullfile(dirName, 's?_dic.tif');     % ORIGINAL CODE
sliceFiles = dir(sliceFilePattern);

for k = 1 : length(sliceFiles)
    baseFileName = sliceFiles(k).name;
    fullFileName = fullfile(sliceFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    
    % store file name and slice number
    sliceNumber = baseFileName(2);
    sliceNumberAsNum = str2num(sliceNumber);
    sliceImgsNames = [sliceImgsNames; fullFileName];  
    
    % read image and store it
    % (I am calling the slice image 'fixed' for registration purposes)
    fixed = imread(fullFileName);
    sliceImgs = [sliceImgs; fixed];
    
    % display slice image in an appropriately named figure
    figure('Name', strcat(expDate, '_', mouseNumber, '_', baseFileName(1:2)));
    imshow(fixed);     
    
    % Now get cell images from each slice
    % "*" is a wildcard/placeholder for one or more characters
    cellFilePattern = fullfile(dirName, strcat('s', sliceNumber, '*x_dic.tif'));
    cellFiles = dir(cellFilePattern);
    
    dataPerSlice = [];
    
    if ~isempty(cellFiles)    
        for i = 1 : length(cellFiles)
            baseFileNameCell = cellFiles(i).name;
            fullFileNameCell = fullfile(cellFiles(i).folder, baseFileNameCell);
            fprintf(1, 'Now reading %s\n', fullFileNameCell);

            % store file name and cell number
            cellNumber = baseFileNameCell(4);
            cellNumberAsNum = str2num(cellNumber);
            cellImgsNames = [cellImgsNames; fullFileNameCell];  

            % read image and store it
            % (I am calling the cell image 'moving' for registration purposes)
            moving = imread(fullFileNameCell);
            cellImgs = [cellImgs; moving];

            % display slice image in an appropriately named figure
            figure('Name', strcat(expDate, '_', mouseNumber, '_', baseFileNameCell(1:4)));
            imshow(moving); 

            % display pipette position (aka cell location) based on user input
            viscircles([pipx pipy], 5, 'Color', 'r', 'LineWidth', 5);

            % collect registration transformation
            tform = imregtform(moving, fixed, 'translation', optimizer, metric);
            xtranslation = tform.T(3,1);
            ytranslation = tform.T(3,2);
            xregistered = pipx + xtranslation;
            yregistered = pipy + ytranslation;

            % align images based on registration transformation
            figure('Name', strcat(expDate, '_', mouseNumber, '_', baseFileNameCell(1:4), " registered to ", baseFileName(1:2)));
            movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));
            imshowpair(fixed, movingRegistered,'Scaling','joint');
            viscircles([xregistered yregistered], 5, 'Color', 'r', 'LineWidth', 5);

            % store all relevant information
            data = [data; expDateAsNum, mouseNumberAsNum, pipx, pipy, sliceNumberAsNum, cellNumberAsNum, xtranslation, ytranslation, xregistered, yregistered];
            dataPerSlice = [dataPerSlice; xregistered, yregistered, cellNumberAsNum];

        end       
                
    % show slice image with all recorded cells
    figure('Name', strcat(expDate, '_', mouseNumber, '_', baseFileName(1:2), '_', 'all recorded cells'));
    imshow(fixed);
    hold on;
    scatter(dataPerSlice(:,1), dataPerSlice(:,2), 600, 'filled', 'MarkerFaceColor', 'white', 'MarkerEdgeColor', 'black', 'LineWidth', 1);
    text(dataPerSlice(:,1)-10, dataPerSlice(:,2)-2, num2str(dataPerSlice(:,3)), 'FontSize', 18, 'Color', 'black');
    hold off; 
    
    end
     
end

% save XLSX file with data 
xlsxFileName = strcat(expDate, '_', mouseNumber, '_', 'cell location');
fulldirectory = strcat(dirName,'\',xlsxFileName,'.xlsx');        
dataInCellFormat = {};
dataInCellFormat = num2cell(data);
labeledData = cell2table(dataInCellFormat,'VariableNames',...
    {'date', 'mouse', 'pipx', 'pipy', 'sliceNumber', 'cellNumber', 'xtranslation', 'ytranslation', 'xregistered', 'yregistered'});
writetable(labeledData,fulldirectory);
disp('I saved it')  

