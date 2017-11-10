function rectMat = rectTextureCreate(width, length, colorIndexRectangle, colorIndexBacground,overlayMat, maskParam, gapLine, noFillIn)
% Create a rectangultar texture containing the matrix defining a rectangle
%
% Input:    width, length : the width and length (pixel)
%           colorIndexRectangle    : the RGB value of the rectangle surface
%           colorIndexBackground              : the RGB value of the background
% Output:   rectMat    : the numerical matrix of the rectangle

% Initiallize the matrix
if isempty(noFillIn) || (noFillIn==0)
    widthMin = 10;
elseif noFillIn
    widthMin = 5;
end
[nRows, nCols, ~] = size(overlayMat);
boundarySize = length;
if ~isempty(maskParam)
    length = length + maskParam + gapLine;
end
rectWidth = max([width nRows+2 widthMin]);
rectLength = max([length+2 nCols widthMin]);
rectMatGray = zeros(rectWidth, rectLength);
rectMatGray(round((rectWidth-width)/2):width+round((rectWidth-width)/2), ...
            round((rectLength-length)/2):length+round((rectLength-length)/2)) = 255;

% Convert to color rectangle
rectMat = NaN(size(rectMatGray,1), size(rectMatGray,2), 3);
rectTemp = NaN(size(rectMatGray,1), size(rectMatGray,2));
for ii = 1 : 3
    rectTemp(rectMatGray == 255) = 255 * colorIndexRectangle(ii);
    rectTemp(rectMatGray == 0) = floor(255 * colorIndexBacground(ii));
    if ~isempty(maskParam)
        maskStart = round(boundarySize/2);
        rowInd = rectTemp(:,10) == 255*colorIndexRectangle(ii);
        rectTemp(rowInd, maskStart+4 : maskStart+maskParam+gapLine-1) = floor(255 * colorIndexBacground(ii));
    end
    rectMat(:,:,ii) = rectTemp;
end

% % Smooth and crop the rectangle
% myfilter = fspecial('gaussian',[2 2], 0.5);
% smoothRect = imfilter(rectMat, myfilter, 'replicate');
% rectMat = smoothRect;

% Overlay another texture
tempMat = rectMat;
if ~isempty(overlayMat)
    xCenter = round(size(tempMat,1)/2);
    xStart = xCenter-round(nRows/2); 
    yCenter = round(size(tempMat,2)/2);
    yStart = yCenter-round(nCols/2);
    tempMat(xStart+1:xStart+nRows, yStart:yStart+nCols-1,:) = overlayMat;
end
tempMat(tempMat>255) = 255;
tempMat(tempMat<0) = 0;
rectMat = round(tempMat);
    
