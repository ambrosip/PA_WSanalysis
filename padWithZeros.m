function paddedMatrix = padWithZeros(centerMatrix, outputSize)
% PADWITHZEROS pads a matrix centerMatrix with zeros all around it. Returns
% the padded matrix that has size given by outputSize.
%
% Syntax
%   paddedMatrix = PADWITHZEROS(centerMatrix, outputSize)

% Vector that stores the x and y paddings size
% padding(1) will give the x padding and padding(2) the y padding
padding = floor((outputSize - size(centerMatrix)) / 2);

% Create a zero matrix of the correct size
paddedMatrix = zeros(outputSize);

% Replaces the center values with the entries in centerMatrix
paddedMatrix(padding(1) + 1: padding(1) + size(centerMatrix, 1), ...
    padding(2) + 1: padding(2) + size(centerMatrix, 2)) = centerMatrix;
end