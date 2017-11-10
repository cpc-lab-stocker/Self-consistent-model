function ellipseMat = ellipseTextureCreate(radiusAxis1,radiusAxis2,angleEllipse,colorIndexEllipse,colorIndexBacground,paramFilter)
% Create an ellipse texture containing the matrix defining an ellipse
%
% Input:    radiusAxis1, radiusAxis2: the radii of 2 axes in pixel
%           angleEllipse            : the ellipse orientation (degree)
%           colorIndexEllipse       : the RGB value of the ellipse
%           colorIndexBackground    : the RGB value of the background
%           paramFilter(3d vector)  : the Gaussian width, height and std
% Output:   ellipseMat              : the numerical matrix of the ellipse
%
if nargin < 6
    paramFilter = [8 8 3];
end

% Create an upsampled binary ellipse structure and then downsample it for
% smoother ellipse
upsampleFactor = 2;
ellipseBinary = Ellipse(upsampleFactor*radiusAxis1,upsampleFactor*radiusAxis2);
ellipseBinary = imresize(ellipseBinary, 1/upsampleFactor, 'bicubic');
ellipseBinary = padarray(ellipseBinary, [10 10], 0);
ellipseGray = ellipseBinary * 255;

% Rotate the ellipse 
if angleEllipse ~= 0
    ellipseGrayRotate = imrotate(ellipseGray, angleEllipse);
else
    ellipseGrayRotate = ellipseGray;
end

% Convert to color ellipse
ellipseMat = NaN(size(ellipseGrayRotate,1), size(ellipseGrayRotate,2), 3);
ellipseTemp = NaN(size(ellipseGrayRotate,1), size(ellipseGrayRotate,2));
for ii = 1 : 3
    ellipseTemp(ellipseGrayRotate == 255) = 255 * colorIndexEllipse(ii);
    ellipseTemp(ellipseGrayRotate == 0) = 255 * colorIndexBacground(ii);    
    ellipseMat(:,:,ii) = ellipseTemp;
end

% Smooth the ellipse
myfilter = fspecial('gaussian', paramFilter(1:2), paramFilter(3));
smoothEllipse = imfilter(ellipseMat, myfilter, 'replicate');
ellipseMat = round(smoothEllipse);

% Tighten the boundary
cutoffVal = 150;
indicatorMat = ellipseMat > cutoffVal;
indicatorMat1 = indicatorMat(:,:,1);
indicatorMat2 = indicatorMat(:,:,2);
indicatorMat3 = indicatorMat(:,:,3);
[~, indexMax] = max([sum(indicatorMat1(:)) sum(indicatorMat2(:)) sum(indicatorMat3(:))]);
indicatorChosen = squeeze(indicatorMat(:,:,indexMax));
xLim = [find(indicatorChosen(round(size(indicatorChosen,1)/2),:), 1, 'first') ...
        find(indicatorChosen(round(size(indicatorChosen,1)/2),:), 1, 'last')];
yLim = [find(indicatorChosen(:,round(size(indicatorChosen,2)/2)), 1, 'first') ...
        find(indicatorChosen(:,round(size(indicatorChosen,2)/2)), 1, 'last')];
ellipseMat = imcrop(ellipseMat,[xLim(1) yLim(1) xLim(2)-xLim(1) yLim(2)-yLim(1)]);    