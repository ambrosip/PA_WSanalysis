saveDir = '/Users/priscilla/OHSU Dropbox/Priscilla Ambrosi/Dropbox - Lerner Lab/Ambrosi et al_sCRACM 2024/Data Analysis/2024-10-28';
databaseFile = "/Users/priscilla/OHSU Dropbox/Priscilla Ambrosi/Dropbox - Lerner Lab/Ambrosi et al_sCRACM 2024/2024-10-28 data for matlab.xlsx";
database = readtable(databaseFile);

rows = height(database);

for row = 1:rows

    cellImageDir = cell2mat(database.dir(row));
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


