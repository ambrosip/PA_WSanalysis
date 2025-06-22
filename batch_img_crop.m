saveDir = '/Users/priscilla/OHSU Dropbox/Priscilla Ambrosi/Dropbox - Lerner Lab/Ambrosi et al_sCRACM 2024/Data Analysis/2024-11-14 renamed cropped imgs';
databaseFile = "/Users/priscilla/OHSU Dropbox/Priscilla Ambrosi/Dropbox - Lerner Lab/Ambrosi et al_sCRACM 2024/2024-10-28 data for matlab.xlsx";
database = readtable(databaseFile);
renameFiles = 1;

rows = height(database);

for row = 1:rows

    cellImageDir = cell2mat(database.dir(row));
    leftCrop = database.leftCrop(row);
    rightCrop = database.rightCrop(row);
    topCrop = database.topCrop(row);
    bottomCrop = database.bottomCrop(row);
    gridColumns = 5;
    gridRows = 5;

    cellNamePrefix = "";
    if renameFiles == 1
        % collect info from the database about this cell
        circuitLabel = cell2mat(database.circuit(row));
        cellLabel = cell2mat(database.cell(row));
        % adjust mouse number according to length to get consistent m000
        % nomenclature
        if length(num2str(database.m(row))) == 3
            mouseNumber = strcat('m',num2str(database.m(row)));
        elseif length(num2str(database.m(row))) == 2
            mouseNumber = strcat('m0',num2str(database.m(row)));
        elseif isscalar(num2str(database.m(row)))
            mouseNumber = strcat('m00',num2str(database.m(row)));  
        end

        % create name prefix we will use to label figures and spreadsheets
        cellNamePrefix = strcat(circuitLabel, '_', mouseNumber, '_', cellLabel, '_');
    end

    cellImageFile = cell2mat(database.img5x5pre(row));
    % commenting out the prefix definition before adding the varible renameFiles
    % prefix = strcat(cell2mat(database.cell(row)), '_5x5pre_');
    prefix = "5x5pre";
    polygon_img_crop(cellImageFile,cellImageDir,leftCrop,rightCrop,topCrop,bottomCrop,gridColumns,gridRows,prefix,renameFiles,cellNamePrefix);

    cellImageFile = cell2mat(database.img5x5post(row));
    % prefix = strcat(cell2mat(database.cell(row)), '_5x5post_');
    prefix = "5x5post";    
    polygon_img_crop(cellImageFile,cellImageDir,leftCrop,rightCrop,topCrop,bottomCrop,gridColumns,gridRows,prefix,renameFiles,cellNamePrefix);

    cellImageFile = cell2mat(database.tracing(row));
    % prefix = strcat(cell2mat(database.cell(row)), '_tracing_');
    prefix = "tracing";
    polygon_img_crop(cellImageFile,cellImageDir,leftCrop,rightCrop,topCrop,bottomCrop,gridColumns,gridRows,prefix,renameFiles,cellNamePrefix);

    if database.img11x7pre(row) ~= "n"
        gridColumns = 11;
        gridRows = 7;
        cellImageFile = cell2mat(database.img11x7pre(row));
        % prefix = strcat(cell2mat(database.cell(row)), '_11x7pre_');
        prefix = "11x7pre";
        polygon_img_crop(cellImageFile,cellImageDir,leftCrop,rightCrop,topCrop,bottomCrop,gridColumns,gridRows,prefix,renameFiles,cellNamePrefix);        
    end

    if database.img11x7post(row) ~= "n"
        gridColumns = 11;
        gridRows = 7;
        cellImageFile = cell2mat(database.img11x7post(row));
        % prefix = strcat(cell2mat(database.cell(row)), '_11x7post_');
        prefix = "11x7post";
        polygon_img_crop(cellImageFile,cellImageDir,leftCrop,rightCrop,topCrop,bottomCrop,gridColumns,gridRows,prefix,renameFiles,cellNamePrefix);
    end

    if database.tracing2(row) ~= "n"
        gridColumns = 11;
        gridRows = 7;
        cellImageFile = cell2mat(database.tracing2(row));
        % prefix = strcat(cell2mat(database.cell(row)), '_tracing2_');
        prefix = "tracing2";
        polygon_img_crop(cellImageFile,cellImageDir,leftCrop,rightCrop,topCrop,bottomCrop,gridColumns,gridRows,prefix,renameFiles,cellNamePrefix);
    end

    saveAllFigs(saveDir);
    close all;

end


