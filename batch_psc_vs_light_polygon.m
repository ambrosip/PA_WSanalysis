% Things you need added to the MATLAB path:
% wavesurfer code
% PA_WSanalysis code

saveDir = '/Users/priscilla/OHSU Dropbox/Priscilla Ambrosi/Dropbox - Lerner Lab/Ambrosi et al_sCRACM 2024/Data Analysis/2024-10-29';
databaseFile = "/Users/priscilla/OHSU Dropbox/Priscilla Ambrosi/Dropbox - Lerner Lab/Ambrosi et al_sCRACM 2024/2024-10-28 data for matlab.xlsx";
database = readtable(databaseFile);

% Polygon designs - sub for orderOfROIs
% this is the order of the design 5x5 spaced out
spaced5x5 = [8 16 14 23 3 10 12 25 21 19 6 18 1 5 11 4 22 15 13 7 2 20 24 17 9]';
% this is the order of the design 9x9 spaced out
spaced9x9 = [50 10 46 73 12 25 69 27 2 4 48 38 40 79 42 67 75 77 44 34 14 32 36 59 71 53 65 6 22 64 16 30 81 57 55 63 8 20 18 28 1 51 9 23 61 78 54 68 76 3 47 49 21 7 24 60 26 66 11 45 39 33 5 62 80 56 70 74 41 43 15 17 52 31 13 58 72 35 37 19 29]';
% this is the order of the design 6x6 max sep
spaced6x6 = [18 31 33 9 35 5 16 22 27 3 7 25 14 29 20 1 11 15 12 24 36 17 13 23 2 4 8 32 28 21 6 34 10 26 19 30]';
% this is the order of the design 11x7 random1
random11x7v1 = [33 65 9 38 46 68 48 62 20 51 47 53 75 35 30 69 71 13 6 34 16 2 44 28 27 73 31 61 39 21 22 8 58 19 63 42 52 66 1 72 26 5 25 18 45 64 29 3 60 56 70 76 41 14 10 32 23 54 15 50 36 67 24 77 37 40 74 57 43 59 17 4 55 12 11 7 49]';

% % Quality check of Polygon designs 
% reshape(spaced5x5,5,5).';
% reshape(spaced9x9,9,9).';
% reshape(spaced6x6,6,6).';
% reshape(random11x7v1,11,7).';

% each row is a cell
rows = height(database);

% iterate through rows (aka cells) and gather info from database
for row = 1:rows

    cellDir = cell2mat(database.dir(row));
    cellObj = WSfileFromDir(cellDir);
    leftCrop = database.leftCrop(row);
    rightCrop = database.rightCrop(row);
    topCrop = database.topCrop(row);
    bottomCrop = database.bottomCrop(row);

    % adjust mouse number according to length to get consistent m000
    % nomenclature
    if length(num2str(database.m(row))) == 3
        mouseNumber = strcat('m',num2str(database.m(row)));
    elseif length(num2str(database.m(row))) == 2
        mouseNumber = strcat('m0',num2str(database.m(row)));
    elseif isscalar(num2str(database.m(row)))
        mouseNumber = strcat('m00',num2str(database.m(row)));  
    end

    cellNamePrefix = strcat(cell2mat(database.circuit(row)), '_', mouseNumber, '_', cell2mat(database.cell(row)), '_');

    sweepPrefix = 'coolLED_';

    sweepPrefix = '5x5_40x_'
    gridColumns = 5;
    gridRows = 5;





% % % % KEEP GOING


    leftCrop = database.leftCrop(row);
    rightCrop = database.rightCrop(row);
    topCrop = database.topCrop(row);
    bottomCrop = database.bottomCrop(row);
    gridColumns = 5;
    gridRows = 5;

    cellImageFile = cell2mat(database.img5x5pre(row));
    prefix = strcat(cell2mat(database.cell(row)), '_5x5pre_');
    polygon_img_crop(cellImageFile,cellImageDir,leftCrop,rightCrop,topCrop,bottomCrop,gridColumns,gridRows,prefix);

    cellImageFile = cell2mat(database.img5x5post(row));
    prefix = strcat(cell2mat(database.cell(row)), '_5x5post_');
    polygon_img_crop(cellImageFile,cellImageDir,leftCrop,rightCrop,topCrop,bottomCrop,gridColumns,gridRows,prefix);

    cellImageFile = cell2mat(database.tracing(row));
    prefix = strcat(cell2mat(database.cell(row)), '_tracing_');
    polygon_img_crop(cellImageFile,cellImageDir,leftCrop,rightCrop,topCrop,bottomCrop,gridColumns,gridRows,prefix);

    if database.img11x7pre(row) ~= "n"
        gridColumns = 11;
        gridRows = 7;
        cellImageFile = cell2mat(database.img11x7pre(row));
        prefix = strcat(cell2mat(database.cell(row)), '_11x7pre_');
        polygon_img_crop(cellImageFile,cellImageDir,leftCrop,rightCrop,topCrop,bottomCrop,gridColumns,gridRows,prefix);        
    end

    if database.img11x7post(row) ~= "n"
        gridColumns = 11;
        gridRows = 7;
        cellImageFile = cell2mat(database.img11x7post(row));
        prefix = strcat(cell2mat(database.cell(row)), '_11x7post_');
        polygon_img_crop(cellImageFile,cellImageDir,leftCrop,rightCrop,topCrop,bottomCrop,gridColumns,gridRows,prefix);
    end

    if database.tracing2(row) ~= "n"
        gridColumns = 11;
        gridRows = 7;
        cellImageFile = cell2mat(database.tracing2(row));
        prefix = strcat(cell2mat(database.cell(row)), '_tracing2_');
        polygon_img_crop(cellImageFile,cellImageDir,leftCrop,rightCrop,topCrop,bottomCrop,gridColumns,gridRows,prefix);
    end

    saveAllFigs(saveDir);
    close all;

end


